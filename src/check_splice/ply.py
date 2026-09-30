import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import py2bit
from Bio.Seq import Seq
from plotly.subplots import make_subplots

from .common import filter_cpcdh_cluster_exon
from .interact import BedpeJustIntronFilter
from .utils import clone2assemble, clone2treat, get_merge_bam, treat2assemble


def get_hover(
    cfg: dict,
    assemble: str,
    query_name: str,
    R1R2: bool,
    query: str,
    ref_blocks: str,
    query_blocks: str,
    align_strings: str,
    idx: int,
) -> str:
    ref_chrom, ref_start, ref_end, ref_strand = ref_blocks.split(";")[idx].split(":")
    ref_start, ref_end = int(ref_start), int(ref_end)
    with py2bit.open(cfg[assemble]["2bit"]) as tb:
        ref_seg = tb.sequence(ref_chrom, ref_start, ref_end)
    if ref_strand == "-":
        ref_seg = str(Seq(ref_seg).reverse_complement())

    query_start, query_end = query_blocks.split(";")[idx].split(":")
    query_start, query_end = int(query_start), int(query_end)
    query_seg = query[query_start:query_end]
    align_string = align_strings.split(";")[idx]

    pattern = re.compile(r"(\d+)([MID])")
    refline = []
    midline = []
    queryline = []
    ref_pos = 0
    query_pos = 0
    for length, op in pattern.findall(align_string):
        length = int(length)
        if op == "I":
            refline.append("-" * length)
            midline.append("." * length)
            queryline.append(query_seg[query_pos : query_pos + length])
            query_pos += length
        elif op == "D":
            refline.append(ref_seg[ref_pos : ref_pos + length])
            midline.append("." * length)
            queryline.append("-" * length)
            ref_pos += length
        else:
            # op == "M"
            refline.append(ref_seg[ref_pos : ref_pos + length])
            midline.append(
                "".join([
                    "|" if ref_seg[ref_pos + i] == query_seg[query_pos + i] else "."
                    for i in range(length)
                ])
            )
            queryline.append(query_seg[query_pos : query_pos + length])
            ref_pos += length
            query_pos += length

    with py2bit.open(cfg[assemble]["2bit"]) as tb:
        if ref_strand == "+":
            refprefix = (
                tb.sequence(ref_chrom, ref_start - query_start, ref_start)
                if query_start > 0
                else ""
            )
            refsuffix = (
                tb.sequence(ref_chrom, ref_end, ref_end + len(query) - query_end)
                if len(query) > query_end
                else ""
            )
        else:
            refprefix = (
                str(
                    Seq(
                        tb.sequence(ref_chrom, ref_end, ref_end + query_start)
                    ).reverse_complement()
                )
                if query_start > 0
                else ""
            )
            refsuffix = (
                str(
                    Seq(
                        tb.sequence(
                            ref_chrom, ref_start - len(query) + query_end, ref_start
                        )
                    ).reverse_complement()
                )
                if len(query) > query_end
                else ""
            )

    refline = refprefix + "".join(refline) + refsuffix
    midline = "." * query_start + "".join(midline) + "." * (len(query) - query_end)
    queryline = query[:query_start] + "".join(queryline) + query[query_end:]

    if ref_strand == "+":
        hover_text = f"{queryline}<br>{midline}<br>{refline}<br>{'|' * len(refline)}<br>{str(Seq(refline).reverse_complement())}<br><br>"
    else:
        hover_text = f"<br><br>{str(Seq(refline).reverse_complement())}<br>{'|' * len(refline)}<br>{refline[::-1]}<br>{midline[::-1]}<br>{queryline[::-1]}"
    hover_text = f"{query_name}<br>{R1R2}<br>{align_string}<br>{hover_text}"

    return hover_text


def get_hovers(
    cfg: dict,
    assemble: str,
    query_name: str,
    R1R2: bool,
    query: str,
    ref_blocks: str,
    query_blocks: str,
    align_strings: str,
) -> str:
    return ";;".join([
        get_hover(
            cfg=cfg,
            assemble=assemble,
            query_name=query_name,
            R1R2=R1R2,
            query=query,
            ref_blocks=ref_blocks,
            query_blocks=query_blocks,
            align_strings=align_strings,
            idx=idx,
        )
        for idx in range(len(ref_blocks.split(";")))
    ])


def get_reads_with_details(cfg: dict) -> pd.DataFrame:
    df = (
        pd
        .read_feather(cfg["data_dir"] / "result" / "reads.feather")
        .query("ref_blocks.str.contains(';')")
        .reset_index(drop=True)
    )

    bedpe_just_intron_filter = BedpeJustIntronFilter(cfg)
    df_expand = (
        pd
        .read_csv(
            cfg["data_dir"] / "result" / "expand_splice.csv",
            header=0,
            keep_default_na=False,
        )
        .assign(
            just=lambda df: bedpe_just_intron_filter(
                starts=np.minimum(df["end1"], df["end2"]),
                ends=np.maximum(df["start1"], df["start2"]),
                assembles=df["clone"].map(clone2assemble),
            ),
            match_percent=lambda df: (
                df["match_base"] / df["joint_block"].str.len() * 100
            ),
            detail=lambda df: df["detail"].str.replace("\n", "<br>"),
        )
        .sort_values(by=["bidx"], ignore_index=True)
        .groupby(
            ["exp", "protein", "clone", "rep", "query_name", "is_read1"], as_index=False
        )
        .agg(
            just=pd.NamedAgg(column="just", aggfunc=any),
            match_percents=pd.NamedAgg(
                column="match_percent", aggfunc=lambda se: ";".join(se.astype(str))
            ),
            details=pd.NamedAgg(column="detail", aggfunc=";;".join),
        )
    )

    df = df.merge(
        right=df_expand[
            [
                "exp",
                "protein",
                "clone",
                "rep",
                "query_name",
                "is_read1",
                "just",
                "match_percents",
                "details",
            ]
        ],
        how="left",
        on=[
            "exp",
            "protein",
            "clone",
            "rep",
            "query_name",
            "is_read1",
        ],
        validate="one_to_one",
    )

    return df


def initialize_fig(
    cfg: dict,
    cluster: str,
    exp: str,
    protein: str,
    treat: str,
    max_ys: pd.Series,
    rows: int,
    block_num: int,
) -> go.Figure:
    max_ys = max_ys.reindex(list(range(1, rows))).fillna(0.0)
    total_bottom_margin = cfg["plotly"]["bottom_margin"] + cfg["plotly"]["range_slider"]
    total_vertical_spacing = (
        cfg["plotly"]["vertical_spacing"] + cfg["plotly"]["range_slider"]
    )
    fig_height = (
        cfg["plotly"]["block_height"] * (block_num + 2 + rows)
        + total_vertical_spacing * (rows - 1)
        + cfg["plotly"]["top_margin"]
        + total_bottom_margin
    )
    row_heights = [bn + 1 for bn in max_ys.to_list() + [2]]
    fig = make_subplots(
        rows=rows,
        cols=1,
        shared_xaxes=True,
        subplot_titles=(
            "just +",
            "just -",
            "strange +",
            "strange -",
            "gene",
        ),
        row_heights=row_heights,
        vertical_spacing=total_vertical_spacing
        / (fig_height - cfg["plotly"]["top_margin"] - total_bottom_margin),
    )
    fig.update_layout(
        title=f"{exp}_{protein}_{treat}",
        height=fig_height,
        margin={
            "t": cfg["plotly"]["top_margin"],
            "b": total_bottom_margin,
        },
        hovermode="closest",
    )
    fig.update_xaxes(
        range=[
            cfg[treat2assemble(treat)][cluster]["start"],
            cfg[treat2assemble(treat)][cluster]["end"],
        ],
        rangeslider={
            "visible": True,
            "thickness": cfg["plotly"]["range_slider"]
            / (
                fig_height
                - cfg["plotly"]["top_margin"]
                - total_bottom_margin
                - total_vertical_spacing * (rows - 1)
            ),
        },
    )
    for row in range(1, 5):
        fig.update_yaxes(range=[-1, max_ys.loc[row]], fixedrange=True, row=row, col=1)
    fig.update_yaxes(range=[-2, 1], row=5, col=1)

    return fig


def get_plotly_interact(cfg: dict, cluster: str) -> None:
    df = (
        get_reads_with_details(cfg)
        .assign(
            strand=lambda df: (
                df["ref_blocks"]
                .str.split(";", n=1, expand=True)[0]
                .str.rsplit(":", n=1, expand=True)[1]
            ),
            row=lambda df: (
                df["just"].map({True: 0, False: 1}) * 2
                + df["strand"].map({"+": 0, "-": 1})
                + 1
            ),
            read_start=lambda df: np.minimum(
                df["ref_blocks"]
                .str.split(";", n=1, expand=True)[0]
                .str.split(":", expand=True)[1]
                .astype(int),
                df["ref_blocks"]
                .str.rsplit(";", n=1, expand=True)[1]
                .str.split(":", expand=True)[1]
                .astype(int),
            ),
            read_end=lambda df: np.maximum(
                df["ref_blocks"]
                .str.split(";", n=1, expand=True)[0]
                .str.split(":", expand=True)[2]
                .astype(int),
                df["ref_blocks"]
                .str.rsplit(";", n=1, expand=True)[1]
                .str.split(":", expand=True)[2]
                .astype(int),
            ),
            read_start_idx=lambda df: list(zip(df["read_start"], df.index)),
            assemble=lambda df: df["clone"].map(clone2assemble),
            treat=lambda df: df["clone"].map(clone2treat),
            cluster_start=lambda df: df["assemble"].map(
                lambda assemble: cfg[assemble][cluster]["start"]
            ),
            cluster_end=lambda df: df["assemble"].map(
                lambda assemble: cfg[assemble][cluster]["end"]
            ),
        )
        .query(
            """
                read_start >= cluster_start and \
                read_end <= cluster_end
            """
        )
        .reset_index(drop=True)
        .assign(
            args=lambda df: (
                df["assemble"]
                + "|"
                + df["query_name"]
                + "|"
                + df["is_read1"].map({True: "R1", False: "R2"})
                + "|"
                + df["query"]
                + "|"
                + df["ref_blocks"]
                + "|"
                + df["query_blocks"]
                + "|"
                + df["align_strings"]
            ),
            hovers=lambda df: df["args"].map(
                lambda args: get_hovers(cfg, *args.split("|"))
            ),
        )
    )

    df_merge = get_merge_bam(cfg)
    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        df_slice = (
            df
            .query("exp == @exp and protein == @protein and treat == @treat")
            .reset_index(drop=True)
            .assign(
                y=lambda df: (
                    df.groupby(["just", "strand"])["read_start_idx"].rank(
                        method="dense"
                    )
                    - 1
                ),
            )
        )

        fig = initialize_fig(
            cfg=cfg,
            cluster=cluster,
            exp=exp,
            protein=protein,
            treat=treat,
            max_ys=df_slice["row"].value_counts(),
            rows=5,
            block_num=len(df_slice),
        )

        for ref_blocks, hovers, match_percents, details, row, y in zip(
            df_slice["ref_blocks"],
            df_slice["hovers"],
            df_slice["match_percents"],
            df_slice["details"],
            df_slice["row"],
            df_slice["y"],
        ):
            ref_blocks = ref_blocks.split(";")
            hovers = hovers.split(";;")
            match_percents = match_percents.split(";")
            details = details.split(";;")
            for bidx, (ref_block, hover) in enumerate(zip(ref_blocks, hovers)):
                ref_chrom, ref_start, ref_end, ref_strand = ref_block.split(":")
                ref_start, ref_end = int(ref_start), int(ref_end)
                fig.add_trace(
                    go.Scatter(
                        x=[ref_start, ref_start, ref_end, ref_end],
                        y=[y + 0.4, y - 0.4, y - 0.4, y + 0.4],
                        mode="none",
                        fill="toself",
                        fillcolor=cfg["plotly"]["fillcolor"]["map"],
                        name=hover,
                        hoverlabel=cfg["plotly"]["hoverlabel"],
                        showlegend=False,
                    ),
                    row=row,
                    col=1,
                )

                if bidx < len(match_percents):
                    next_ref_block = ref_blocks[bidx + 1]
                    match_percent = match_percents[bidx]
                    detail = details[bidx]
                    next_ref_chrom, next_ref_start, next_ref_end, next_ref_strand = (
                        next_ref_block.split(":")
                    )
                    next_ref_start, next_ref_end = (
                        int(next_ref_start),
                        int(next_ref_end),
                    )
                    link_start = min(ref_end, next_ref_end)
                    link_end = max(ref_start, next_ref_start)
                    fig.add_trace(
                        go.Scatter(
                            x=[link_start, link_start, link_end, link_end],
                            y=[y + 0.2, y - 0.2, y - 0.2, y + 0.2],
                            mode="none",
                            fill="toself",
                            fillcolor=cfg["plotly"]["fillcolor"]["splice"],
                            name=f"match percent: {match_percent}<br>{detail}",
                            hoverlabel=cfg["plotly"]["hoverlabel"],
                            showlegend=False,
                        ),
                        row=row,
                        col=1,
                    )

        assemble = treat2assemble(treat)
        df_cpcdh = filter_cpcdh_cluster_exon(
            df_cpcdh=pd.read_csv(
                cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0
            ),
            cluster=cluster,
        )
        for start, end, name in zip(
            df_cpcdh["start"], df_cpcdh["end"], df_cpcdh["name"]
        ):
            fig.add_trace(
                go.Scatter(
                    x=[start, start, end, end],
                    y=[0.5, -0.5, -0.5, 0.5],
                    mode="none",
                    fill="toself",
                    fillcolor=cfg["plotly"]["fillcolor"]["exon"],
                    name=name,
                    hoverlabel=cfg["plotly"]["hoverlabel"],
                    showlegend=False,
                ),
                row=5,
                col=1,
            )
            fig.add_annotation(
                x=(start + end) / 2,
                y=-1,
                text=name,
                font=cfg["plotly"]["font"],
                showarrow=False,
                row=5,
                col=1,
            )

        fig.write_html(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "plotly"
            / f"{exp}_{protein}_{treat}.html"
        )

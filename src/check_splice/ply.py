import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import py2bit
from Bio.Seq import Seq
from plotly.subplots import make_subplots

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


def get_plotly_interact(cfg: dict, cluster: str) -> None:
    df = (
        pd
        .read_feather(cfg["data_dir"] / "result" / "reads.feather")
        .query("ref_blocks.str.contains(';')")
        .reset_index(drop=True)
        .assign(
            assemble=lambda df: df["clone"].map(clone2assemble),
            treat=lambda df: df["clone"].map(clone2treat),
            cluster_chrom=lambda df: df["assemble"].map(
                lambda assemble: cfg[assemble][cluster]["chrom"]
            ),
            cluster_start=lambda df: df["assemble"].map(
                lambda assemble: cfg[assemble][cluster]["start"]
            ),
            cluster_end=lambda df: df["assemble"].map(
                lambda assemble: cfg[assemble][cluster]["end"]
            ),
        )
        .query(
            """
                chrom1 == cluster_chrom and \
                start1 >= cluster_start and \
                end1 <= cluster_end and \
                chrom2 == cluster_chrom and \
                start2 >= cluster_start and \
                end2 <= cluster_end
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

    bedpe_just_intron_filter = BedpeJustIntronFilter(cfg)
    df_expand = (
        pd
        .read_csv(cfg["data_dir"] / "result" / "expand_splice.csv", header=0)
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
        .sort_values(
            by=["exp", "protein", "clone", "rep", "query_name", "is_read1", "bidx"]
        )
        .groupby(
            ["exp", "protein", "clone", "rep", "query_name", "is_read1"], ax_index=False
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
            "exp",
            "protein",
            "clone",
            "rep",
            "query_name",
            "is_read1",
            "just",
            "match_percents",
            "details",
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

    df_merge = get_merge_bam(cfg)
    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        df_slice = df.query(
            "exp == @exp and protein == @protein and treat == @treat"
        ).reset_index(drop=True)

        idxs = []
        ref_starts = []
        ref_ends = []
        ref_strands = []
        sam_texts = []
        justs = []
        blat_prevs = []
        blat_nexts = []
        for idx, (
            ref_blocks,
            hovers,
            just,
            match_percents,
            details,
        ) in enumerate(
            zip(
                df_slice["ref_blocks"],
                df_slice["hovers"],
                df_slice["just"],
                df_slice["match_percents"],
                df_slice["details"],
            )
        ):
            for ref_block, hover in zip(ref_blocks.split(";"), hovers.split(";;")):
                ref_chrom, ref_start, ref_end, ref_strand = ref_block.split(":")
                ref_start, ref_end = int(ref_start), int(ref_end)

                idxs.append(idx)
                ref_starts.append(ref_start)
                ref_ends.append(ref_end)
                ref_strands.append(ref_strand)
                sam_texts.append(hover)
                justs.append(just)

            blat_prevs.append("N/A")
            for match_percent, detail in zip(
                match_percents.split(";"), details.split(";;")
            ):
                blat_string = f"match_percent: {match_percent}<br>{detail}"
                blat_prevs.append(blat_string)
                blat_nexts.append(blat_string)
            blat_nexts.append("N/A")

        df_format = (
            pd
            .DataFrame({
                "idx": idxs,
                "ref_start": ref_starts,
                "ref_end": ref_ends,
                "ref_strand": ref_strands,
                "sam_text": sam_texts,
                "just": justs,
                "blat_prev": blat_prevs,
                "blat_next": blat_nexts,
            })
            .sort_values(by=["ref_start", "ref_end"], ignore_index=True)
            .assign(
                just_idx=lambda df: df["just"].map({True: 0, False: 1}),
                ref_strand_idx=lambda df: df["ref_strand"].map({"-": 0, "+": 1}),
                y=lambda df: df.groupby([
                    "just_idx",
                    "ref_strand_idx",
                    "idx",
                ]).transform("ngroup"),
            )
        )

        fig = make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            subplot_titles=(
                "splice",
                "gene",
            ),
            row_heights=[0.95, 0.05],
        )

        for (
            ref_start,
            ref_end,
            y,
            ref_strand,
            sam_text,
            just,
            blat_prev,
            blat_next,
        ) in zip(
            df_format["ref_start"],
            df_format["ref_end"],
            df_format["y"],
            df_format["ref_strand"],
            df_format["sam_text"],
            df_format["just"],
            df_format["blat_prev"],
            df_format["blat_next"],
        ):
            hovertext = f"""
<div>
{sam_text}
</div>
<div style="{{display: flex; gap: 20px;}}">
    <div style="{{flex: 1}}">
        {blat_prev}
    </div>
    <div style="{{flex: 1}}">
        {blat_next}
    </div>
</div>
            """
            fig.add_trace(
                go.Scatter(
                    x=[ref_start, ref_start, ref_end, ref_end],
                    y=[y + 0.5, y - 0.5, y - 0.5, y + 0.5],
                    fill="toself",
                    fillcolor=cfg["color"][ref_strand][just],
                    mode="none",
                    name=hovertext,
                    hoverlabel={
                        "font": {
                            "family": "Courier New, monospace",
                            "size": 14,
                            "color": "black",
                        },
                        "bgcolor": "white",
                    },
                    showlegend=False,
                ),
                row=1,
                col=1,
            )

        assemble = treat2assemble(treat)
        df_cpcdh = pd.read_csv(
            cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0
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
                    fillcolor="blue",
                    name=name,
                    hoverlabel={
                        "font": {
                            "family": "Courier New, monospace",
                            "size": 14,
                            "color": "black",
                        },
                        "bgcolor": "white",
                    },
                    showlegend=False,
                ),
                row=2,
                col=1,
            )
            fig.add_annotation(
                x=(start + end) / 2,
                y=-1,
                text=name,
                font={
                    "size": 8,
                },
                showarrow=False,
                row=2,
                col=1,
            )

        fig.update_layout(
            title=f"{exp}_{protein}_{treat}",
            hovermode="closest",
        )
        fig.update_xaxes(
            range=[
                cfg[treat2assemble(treat)][cluster]["start"],
                cfg[treat2assemble(treat)][cluster]["end"],
            ]
        )
        fig.update_yaxes(range=[-1, len(df_slice)], row=1, col=1)
        fig.update_yaxes(range=[-2, 1], row=2, col=1)

        fig.write_html(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "plotly"
            / f"{exp}_{protein}_{treat}.html"
        )

import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import py2bit
from Bio.Seq import Seq
from plotly.subplots import make_subplots

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
    if query_name == "A00358:943:H7FKWDSX7:2:2313:9787:16939" and R1R2 == "R1":
        breakpoint()
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


def in_range(cfg: dict, cluster: str, assemble: str, ref_blocks: str) -> bool:
    chrom = cfg[assemble][cluster]["chrom"]
    start = cfg[assemble][cluster]["start"]
    end = cfg[assemble][cluster]["end"]
    for ref_block in ref_blocks.split(";"):
        ref_chrom, ref_start, ref_end, ref_strand = ref_block.split(":")
        ref_start, ref_end = int(ref_start), int(ref_end)
        if ref_chrom == chrom and ref_start >= start and ref_end <= end:
            return True

    return False


def get_plotly_interact(cfg: dict, cluster: str) -> None:
    df = (
        (
            pd
            .read_feather(cfg["data_dir"] / "result" / "reads.feather")
            .query("ref_blocks.str.contains(';')")
            .reset_index(drop=True)
        )
        .assign(
            treat=lambda df: df["clone"].map(clone2treat),
            assemble=lambda df: df["clone"].map(clone2assemble),
            assemble_ref_blocks=lambda df: df["assemble"] + "|" + df["ref_blocks"],
            in_cluster=lambda df: df["assemble_ref_blocks"].map(
                lambda assemble_ref_blocks: in_range(
                    cfg, cluster, *assemble_ref_blocks.split("|")
                )
            ),
        )
        .query("in_cluster")
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
        df_slice = df.query(
            "exp == @exp and protein == @protein and treat == @treat"
        ).reset_index(drop=True)
        fig = go.Figure()

        for y, (ref_blocks, hovers) in enumerate(
            zip(df_slice["ref_blocks"], df_slice["hovers"])
        ):
            for ref_block, hover in zip(ref_blocks.split(";"), hovers.split(";;")):
                ref_chrom, ref_start, ref_end, ref_strand = ref_block.split(":")
                ref_start, ref_end = int(ref_start), int(ref_end)

                fig.add_trace(
                    go.Scatter(
                        x=[
                            ref_start,
                            ref_start,
                            ref_end,
                            ref_end,
                            ref_start,
                            (ref_start + ref_end) / 2,
                        ],
                        y=[y - 0.5, y + 0.5, y + 0.5, y - 0.5, y - 0.5, y],
                        fill="toself",
                        mode="text",
                        text=["", "", "", "", "", "+"],
                        hovertext=hover,
                        hovertemplate="%{text}",
                        hoverlabel={
                            "font": {
                                "family": "Courier New, monospace",
                                "size": 14,
                            },
                            "bgcolor": "rgba(255, 255, 255, 0.0)",
                        },
                        showlegend=False,
                    )
                )

        fig.update_layout(
            title=f"{exp}_{protein}_{treat}",
            xaxis={
                "range": [
                    cfg[treat2assemble(treat)][cluster]["start"],
                    cfg[treat2assemble(treat)][cluster]["end"],
                ]
            },
            yaxis={"range": [-1, y + 1]},
            hovermode="closest",
        )

        fig.write_html(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "plotly"
            / f"{exp}_{protein}_{treat}.html"
        )

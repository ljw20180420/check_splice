import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import py2bit
import pysam
from Bio.Seq import Seq

from .utils import get_treat, select_total_count


def get_tas(ref_block: str, query_block: str, align_string: str):
    pattern = re.compile(r"(\d+)([MID])")
    ref = []
    mid = []
    query = []
    ref_pos = 0
    query_pos = 0
    for length, op in pattern.findall(align_string):
        length = int(length)
        if op == "I":
            ref.append("-" * length)
            mid.append("*" * length)
            query.append(query_block[query_pos : query_pos + length])
            query_pos += length
        elif op == "D":
            ref.append(ref_block[ref_pos : ref_pos + length])
            mid.append("*" * length)
            query.append("-" * length)
            ref_pos += length
        else:
            # op == "M"
            ref.append(ref_block[ref_pos : ref_pos + length])
            mid.append(
                "".join([
                    "|" if ref_block[ref_pos + i] == query_block[query_pos + i] else " "
                    for i in range(length)
                ])
            )
            query.append(query_block[query_pos : query_pos + length])
            ref_pos += length
            query_pos += length

    return ref, mid, query


def get_plotly_interact(cfg: dict, assemble: str) -> None:
    def adjacent_block(row: str) -> list:
        blocks = row["blocks"].split(";")
        query_blocks = row["query_blocks"].split(";")
        align_strings = row["align_strings"].split(";")
        interact = []
        for i in range(len(blocks) - 1):
            chrom, ref_start, ref_end, strand = blocks[i].split(":")
            next_chrom, next_ref_start, next_ref_end, next_strand = blocks[i + 1].split(
                ":"
            )
            query_start, query_end = query_blocks[i].split(":")
            next_query_start, next_query_end = query_blocks[i + 1].split(":")
            align_string = align_strings[i]
            next_align_string = align_strings[i + 1]
            assert chrom == next_chrom, "chrom is not consistent"

            if ref_start > next_ref_start:
                chrom, next_chrom = next_chrom, chrom
                ref_start, next_ref_start = next_ref_start, ref_start
                ref_end, next_ref_end = next_ref_end, ref_end
                strand, next_strand = next_strand, strand
                query_start, next_query_start = next_query_start, query_start
                query_end, next_query_end = next_query_end, query_end
                align_string, next_align_string = next_align_string, align_string

            interact.append(
                f"{chrom}:{ref_start}:{ref_end}:{query_start}:{query_end}:{align_string}:{strand}:{next_chrom}:{next_ref_start}:{next_ref_end}:{next_query_start}:{next_query_end}:{next_align_string}:{next_strand}"
            )

        return interact

    def triple_align_string(
        row: pd.Series, ref: str, ref_offset: int, context_size: int
    ):
        def get_ref_prefix_block_suffix(
            ref: str,
            ref_offset: int,
            context_size: int,
            ref_start: int,
            ref_end: int,
            strand: str,
        ) -> tuple[str, str, str]:
            ref_block = ref[ref_start - ref_offset : ref_end - ref_offset]
            ref_prefix = ref[
                ref_start - ref_offset - context_size : ref_start - ref_offset
            ]
            ref_suffix = ref[ref_end - ref_offset : ref_end - ref_offset + context_size]

            if strand == "-":
                ref_block = str(Seq(ref_block).reverse_complement())
                ref_prefix, ref_suffix = (
                    str(Seq(ref_suffix).reverse_complement()),
                    str(Seq(ref_prefix).reverse_complement()),
                )

            return ref_prefix, ref_block, ref_suffix

        def pad_short(
            ref_prefix: str, ref_suffix: str, query_prefix: str, query_suffix: str
        ):
            if len(query_prefix) < len(ref_prefix):
                query_prefix = (
                    "-" * (len(ref_prefix) - len(query_prefix)) + query_prefix
                )
            elif len(query_prefix) > len(ref_prefix):
                ref_prefix = "-" * (len(query_prefix) - len(ref_prefix)) + ref_prefix

            if len(query_suffix) < len(ref_suffix):
                query_suffix = query_suffix + "-" * (
                    len(ref_suffix) - len(query_suffix)
                )
            elif len(query_suffix) > len(ref_suffix):
                ref_suffix = ref_suffix + "-" * (len(query_suffix) - len(ref_suffix))

            mid_prefix = "*" * len(ref_prefix)
            mid_suffix = "*" * len(ref_suffix)

            return (
                ref_prefix,
                ref_suffix,
                mid_prefix,
                mid_suffix,
                query_prefix,
                query_suffix,
            )

        ref_prefix1, ref_block1, ref_suffix1 = get_ref_prefix_block_suffix(
            ref,
            ref_offset,
            context_size,
            row["ref_start1"],
            row["ref_end1"],
            row["strand1"],
        )
        ref_prefix2, ref_block2, ref_suffix2 = get_ref_prefix_block_suffix(
            ref,
            ref_offset,
            context_size,
            row["ref_start2"],
            row["ref_end2"],
            row["strand2"],
        )

        query_block1 = row["query"][row["query_start1"] : row["query_end1"]]
        query_prefix1 = row["query"][: row["query_start1"]]
        query_suffix1 = row["query"][row["query_end1"] :]
        query_block2 = row["query"][row["query_start2"] : row["query_end2"]]
        query_prefix2 = row["query"][: row["query_start2"]]
        query_suffix2 = row["query"][row["query_end2"] :]

        (
            ref_prefix1,
            ref_suffix1,
            mid_prefix1,
            mid_suffix1,
            query_prefix1,
            query_suffix1,
        ) = pad_short(ref_prefix1, ref_suffix1, query_prefix1, query_suffix1)
        (
            ref_prefix2,
            ref_suffix2,
            mid_prefix2,
            mid_suffix2,
            query_prefix2,
            query_suffix2,
        ) = pad_short(ref_prefix2, ref_suffix2, query_prefix2, query_suffix2)

        ref1, mid1, query1 = get_tas(ref_block1, query_block1, row["align_string1"])
        ref2, mid2, query2 = get_tas(ref_block2, query_block2, row["align_string2"])

        refline1 = f"{ref_prefix1}{ref1}{ref_suffix1}"
        midline1 = f"{mid_prefix1}{mid1}{mid_suffix1}"
        queryline1 = f"{query_prefix1}{query1}{query_suffix1}"
        refline2 = f"{ref_prefix2}{ref2}{ref_suffix2}"
        midline2 = f"{mid_prefix2}{mid2}{mid_suffix2}"
        queryline2 = f"{query_prefix2}{query2}{query_suffix2}"

        return pd.Series(
            data=[
                f"{refline1};{midline1};{queryline1}",
                f"{refline2};{midline2};{queryline2}",
            ],
            index=["tas1", "tas2"],
        )

    result_file = cfg["data_dir"] / "result" / f"{assemble}_reads.feather"
    df = pd.read_feather(result_file).assign(
        treat=lambda df: get_treat(df),
        exp_protein_treat=lambda df: (
            df["exp"] + "_" + df["protein"] + "_" + df["treat"]
        ),
    )

    df = (
        df
        .query("not is_shadow and blocks.str.contains(';')")
        .reset_index(drop=True)[
            [
                "exp_protein_treat",
                "query_name",
                "is_read1",
                "blocks",
                "query_blocks",
                "align_strings",
            ]
        ]
        .assign(interact=lambda df: df.apply(adjacent_block, axis=1))
        .explode("interact", ignore_index=True)
    )

    df = (
        pd
        .concat(
            [
                df[["exp_protein_treat", "query_name", "is_read1"]],
                df["interact"]
                .str.split(":", expand=True)
                .rename(
                    columns={
                        0: "chrom1",
                        1: "ref_start1",
                        2: "ref_end1",
                        3: "query_start1",
                        4: "query_end1",
                        5: "align_string1",
                        6: "strand1",
                        7: "chrom2",
                        8: "ref_start2",
                        9: "ref_end2",
                        10: "query_start2",
                        11: "query_end2",
                        12: "align_string2",
                        13: "strand2",
                    }
                )
                .astype({
                    "ref_start1": "int64",
                    "ref_end1": "int64",
                    "query_start1": "int64",
                    "query_end1": "int64",
                    "ref_start2": "int64",
                    "ref_end2": "int64",
                    "query_start2": "int64",
                    "query_end2": "int64",
                }),
            ],
            axis=1,
        )
        .query(
            "ref_start1 >= @cfg[@assemble]['start'] and ref_end2 <= @cfg[@assemble]['end']"
        )
        .sort_values(by=["chrom1", "ref_start1"], ignore_index=True)
    )

    with py2bit.open(cfg[assemble]["2bit"]) as tb:
        ref = tb.sequence(
            cfg[assemble]["chrom"], cfg[assemble]["start"], cfg[assemble]["end"]
        )
        ref_offset = cfg[assemble]["start"]

    df = df.assign(
        motif1=lambda df, ref=ref, ref_offset=ref_offset: df["ref_end1"].map(
            lambda ref_end1: ref[ref_end1 - ref_offset : ref_end1 - ref_offset + 2]
        ),
        motif2=lambda df, ref=ref, ref_offset=ref_offset: df["ref_start2"].map(
            lambda ref_start2: ref[
                ref_start2 - ref_offset - 2 : ref_start2 - ref_offset
            ]
        ),
    )

    df["query"] = float("nan")
    for exp_protein_treat in df["exp_protein_treat"].unique():
        print(exp_protein_treat)
        splice_bam = (
            cfg["data_dir"] / "bam" / "merge" / "splice" / f"{exp_protein_treat}.bam"
        )
        with pysam.AlignmentFile(splice_bam, "rb") as bam:
            bam_index = pysam.IndexedReads(bam)
            bam_index.build()

            df_slice = df.query("exp_protein_treat == @exp_protein_treat")
            indices = []
            queries = []
            for index, row in df_slice.iterrows():
                find_flag = False
                for read in bam_index.find(row["query_name"]):
                    if read.is_unmapped:
                        continue
                    if read.is_secondary:
                        continue
                    if read.is_supplementary:
                        continue
                    if read.is_read1 != row["is_read1"]:
                        continue

                    find_flag = True
                    indices.append(index)
                    queries.append(read.get_forward_sequence())
                    break

                assert find_flag, "not find query name"

        df["query"] = pd.Series(data=queries, index=indices).combine_first(df["query"])

    df = pd.concat(
        (
            df,
            df.apply(
                lambda row, ref=ref, ref_offset=ref_offset, context_size=20: (
                    triple_align_string(row, ref, ref_offset, context_size)
                ),
                axis=1,
            ),
        ),
        axis=1,
    )

    df = df.assign(
        info=lambda df: (
            df["query_name"]
            + df["is_read1"].map({True: "@R1", False: "@R2"})
            + ";"
            + df["motif1"]
            + ":"
            + df["motif2"]
            + ";ref1:"
            + df["ref_start1"]
            + "-"
            + df["ref_end1"]
            + ":"
            + df["strand1"]
            + ";query1:"
            + df["query_start1"]
            + "-"
            + df["query_end1"]
            + ";"
            + df["tas1"]
            + ";ref2:"
            + df["ref_start2"]
            + "-"
            + df["ref_end2"]
            + ":"
            + df["strand2"]
            + ";query2:"
            + df["query_start2"]
            + "-"
            + df["query_end2"]
            + ";"
            + df["tas2"]
        )
    )

    df = (
        df
        .groupby(
            ["exp_protein_treat", "ref_end1", "ref_start2", "motif1", "motif2"],
            as_index=False,
        )
        .agg(
            count=pd.NamedAgg("ref_end1", "size"),
            info=pd.NamedAgg("info", ";;".join),
        )
        .assign(
            total_count=lambda df: df["exp_protein_treat"].map(
                lambda exp_protein_treat, cfg=cfg: select_total_count(
                    cfg, *exp_protein_treat.split("_")
                )
            )
        )
    )

    (cfg["data_dir"] / "result" / "hic" / "plotly").mkdir(exist_ok=True, parents=True)
    df.to_csv(
        cfg["data_dir"] / "result" / "hic" / "plotly" / f"{assemble}_ply.bed",
        sep="\t",
        index=False,
    )


def draw_interact(cfg: dict, assemble: str):
    def create_arc(
        ref_end1: int,
        ref_start2: int,
        num_points=100,
    ):
        angles = np.linspace(0, np.pi, num_points)
        center_x = (ref_end1 + ref_start2) / 2
        radius = (ref_start2 - ref_end1) / 2
        x = center_x + radius * np.cos(angles)
        y = radius * np.sin(angles)
        return x, y, radius

    df = pd.read_csv(
        cfg["data_dir"] / "result" / "hic" / "plotly" / f"{assemble}_ply.bed",
        sep="\t",
        header=0,
    ).assign(
        RPM=lambda df: df["count"] / df["total_count"] * 1_000_000,
    )

    alpha_start = cfg[assemble]["start"]
    alpha_end = cfg[assemble]["end"]

    (cfg["data_dir"] / "result" / "hic" / "plotly" / "draw").mkdir(
        exist_ok=True, parents=True
    )
    for exp_protein_treat in df["exp_protein_treat"].unique():
        df_slice = df.query("exp_protein_treat == @exp_protein_treat").reset_index()
        max_radius = 0
        fig = go.Figure()
        for ref_end1, ref_start2, info, motif1, motif2, RPM in zip(
            df_slice["ref_end1"],
            df_slice["ref_start2"],
            df_slice["info"],
            df_slice["motif1"],
            df_slice["motif2"],
            df_slice["RPM"],
        ):
            if ref_end1 < alpha_start or ref_start2 > alpha_end:
                continue

            if not (motif1 == "GT" and motif2 == "AG") and not (
                motif1 == "CT" and motif2 == "AC"
            ):
                continue

            strand = "+" if motif1 == "GT" and motif2 == "AG" else "-"

            x_coords, y_coords, radius = create_arc(ref_end1, ref_start2)
            if strand == "-":
                y_coords = -y_coords

            max_radius = max(max_radius, radius)

            arc_name = f"{ref_end1}-{ref_start2}:{motif1}-{motif2}"
            fig.add_trace(
                go.Scatter(
                    x=x_coords,
                    y=y_coords,
                    mode="lines",
                    name=arc_name,
                    hovertext=[info.replace(";", "<br>")] * len(x_coords),
                    hovertemplate="%{hovertext}",
                    hoverlabel={
                        "font": {
                            "family": "Courier New, monospace",
                            "size": 14,
                        }
                    },
                    line={
                        "color": "RoyalBlue",
                        "width": 0.5 + RPM,
                    },
                )
            )

        fig.update_layout(
            title=exp_protein_treat,
            xaxis={"range": [alpha_start, alpha_end]},
            yaxis={"range": [-1.1 * max_radius, 1.1 * max_radius]},
            yaxis_scaleanchor="x",
            hovermode="closest",
        )

        fig.write_html(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "plotly"
            / "draw"
            / f"{exp_protein_treat}.html"
        )

import os
import shutil

import numpy as np
import pandas as pd
import py2bit

from .common import (
    BlatSplice,
    donor_acceptor_to_strand,
    get_cpcdh_intron,
    interact2pairs,
    pairs2bedpe,
    substract_bedpe,
    summation_bedpe,
)
from .utils import (
    SelectTotalCount,
    clone2assemble,
    clone2treat,
    get_merge_bam,
    map_to_wild_type_merge,
    treat2assemble,
    treat2diff,
)


def get_interact(cfg: dict) -> None:
    def parse_blocks(ref_query_blocks: str) -> list:
        ref_blocks, query_blocks, query = ref_query_blocks.split("|")
        ref_blocks = ref_blocks.split(";")
        query_blocks = query_blocks.split(";")
        interact = []
        for i in range(len(ref_blocks) - 1):
            ref_block_chrom, ref_block_start, ref_block_end, ref_block_strand = (
                ref_blocks[i].split(":")
            )
            ref_block_start, ref_block_end = int(ref_block_start), int(ref_block_end)
            (
                next_ref_block_chrom,
                next_ref_block_start,
                next_ref_block_end,
                next_ref_block_strand,
            ) = ref_blocks[i + 1].split(":")
            next_ref_block_start, next_ref_block_end = (
                int(next_ref_block_start),
                int(next_ref_block_end),
            )
            assert ref_block_chrom == next_ref_block_chrom, "chrom is not consistent"

            query_block_start, query_block_end = query_blocks[i].split(":")
            query_block_start, query_block_end = (
                int(query_block_start),
                int(query_block_end),
            )
            next_query_block_start, next_query_block_end = query_blocks[i + 1].split(
                ":"
            )
            next_query_block_start, next_query_block_end = (
                int(next_query_block_start),
                int(next_query_block_end),
            )
            joint_block = query[query_block_start:next_query_block_end]

            chrom = ref_block_chrom
            chromStart = min(ref_block_start, next_ref_block_start)
            chromEnd = max(ref_block_end, next_ref_block_end)
            sourceChrom = ref_block_chrom
            sourceStart = ref_block_start
            sourceEnd = ref_block_end
            sourceStrand = ref_block_strand
            targetChrom = next_ref_block_chrom
            targetStart = next_ref_block_start
            targetEnd = next_ref_block_end
            targetStrand = next_ref_block_strand

            interact.append(
                f"{chrom}:{chromStart}:{chromEnd}:{sourceChrom}:{sourceStart}:{sourceEnd}:{sourceStrand}:{targetChrom}:{targetStart}:{targetEnd}:{targetStrand}:{joint_block}:{i}"
            )

        return interact

    df = (
        pd
        .read_feather(cfg["data_dir"] / "result" / "reads.feather")
        .query("ref_blocks.str.contains(';')")
        .reset_index(drop=True)
        .assign(
            treat=lambda df: df["clone"].map(clone2treat),
            assemble=lambda df: df["clone"].map(clone2assemble),
            exp_protein_treat=lambda df: (
                df["exp"] + "_" + df["protein"] + "_" + df["treat"]
            ),
            ref_query_blocks=lambda df: (
                df["ref_blocks"] + "|" + df["query_blocks"] + "|" + df["query"]
            ),
            uid=lambda df: (
                df["assemble"]
                + "_"
                + df["exp"]
                + "_"
                + df["protein"]
                + "_"
                + df["clone"]
                + "_"
                + df["rep"]
                + "_"
                + df["query_name"]
                + "_"
                + df["is_read1"].map({True: "R1", False: "R2"})
            ),
            parse=lambda df: df["ref_query_blocks"].map(parse_blocks),
        )[["assemble", "exp_protein_treat", "query_name", "uid", "parse"]]
        .explode("parse", ignore_index=True)
    )

    df = pd.concat(
        [
            df[["assemble", "exp_protein_treat", "query_name", "uid"]].rename(
                columns={"exp_protein_treat": "exp", "query_name": "name"}
            ),
            df["parse"]
            .str.split(":", expand=True)
            .rename(
                columns={
                    0: "chrom",
                    1: "chromStart",
                    2: "chromEnd",
                    3: "sourceChrom",
                    4: "sourceStart",
                    5: "sourceEnd",
                    6: "sourceStrand",
                    7: "targetChrom",
                    8: "targetStart",
                    9: "targetEnd",
                    10: "targetStrand",
                    11: "joint_block",
                    12: "bidx",
                }
            )
            .astype({
                "chromStart": int,
                "chromEnd": int,
                "sourceStart": int,
                "sourceEnd": int,
                "targetStart": int,
                "targetEnd": int,
            }),
        ],
        axis=1,
    ).assign(uid=lambda df: df["uid"] + "_" + df["bidx"])

    blat_splice = BlatSplice(cfg)
    df_pidents = []
    for assemble in df["assemble"].unique():
        df_slice = df.query("assemble == @assemble")[
            ["uid", "joint_block"]
        ].reset_index(drop=True)
        blat_input = "\n".join(">" + df_slice["uid"] + "\n" + df_slice["joint_block"])
        blat_result = blat_splice(blat_input, assemble)
        df_pident = blat_result.groupby("qseqid", as_index=False)["pident"].max()
        df_pidents.append(df_pident)

    df = df.merge(
        pd.concat(df_pidents, ignore_index=True),
        how="left",
        left_on="uid",
        right_on="qseqid",
        validate="one_to_one",
    ).assign(
        pident=lambda df: df["pident"].fillna(0.0),
    )

    df = df.assign(
        color=0,
        value=1.0,
        sourceName=".",
        targetName=".",
        score=lambda df: 1000 - 10 * df["pident"],
    )[
        [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "value",
            "exp",
            "color",
            "sourceChrom",
            "sourceStart",
            "sourceEnd",
            "sourceName",
            "sourceStrand",
            "targetChrom",
            "targetStart",
            "targetEnd",
            "targetName",
            "targetStrand",
            "assemble",
        ]
    ].sort_values(by=["chrom", "chromStart", "chromEnd"], ignore_index=True)

    df = df.assign(
        minEnd=lambda df: np.minimum(df["sourceEnd"], df["targetEnd"]),
        maxStart=lambda df: np.maximum(df["sourceStart"], df["targetStart"]),
    )

    tbs = {}
    with (
        py2bit.open(cfg["hg19"]["2bit"]) as tbs["hg19"],
        py2bit.open(cfg["mm10"]["2bit"]) as tbs["mm10"],
    ):
        df = df.assign(
            donor=lambda df: [
                tbs[assemble].sequence(chrom, minEnd, minEnd + 2)
                for chrom, minEnd, assemble in zip(
                    df["chrom"], df["minEnd"], df["assemble"]
                )
            ],
            acceptor=lambda df: [
                tbs[assemble].sequence(chrom, maxStart - 2, maxStart)
                for chrom, maxStart, assemble in zip(
                    df["chrom"], df["maxStart"], df["assemble"]
                )
            ],
        )

    df = df.assign(
        sourceStrand=lambda df: donor_acceptor_to_strand(df["donor"], df["acceptor"]),
        targetStrand=lambda df: df["sourceStrand"],
    )[
        [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "value",
            "exp",
            "color",
            "sourceChrom",
            "sourceStart",
            "sourceEnd",
            "sourceName",
            "sourceStrand",
            "targetChrom",
            "targetStart",
            "targetEnd",
            "targetName",
            "targetStrand",
        ]
    ]

    (cfg["data_dir"] / "result" / "hic" / "interact").mkdir(exist_ok=True, parents=True)
    for exp_protein_treat in df["exp"].unique():
        df.query("exp == @exp_protein_treat").to_csv(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "interact"
            / f"{exp_protein_treat}.bed",
            sep="\t",
            header=False,
            index=False,
        )


def interact_to_pairs(cfg: dict) -> None:
    (cfg["data_dir"] / "result" / "hic" / "pairs").mkdir(parents=True, exist_ok=True)
    df_merge = get_merge_bam(cfg)
    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        interact2pairs(
            interact_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "interact"
            / f"{exp}_{protein}_{treat}.bed",
            pairs_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "pairs"
            / f"{exp}_{protein}_{treat}.pairs",
        )


def pairs_to_bedpe(cfg: dict) -> None:
    shutil.rmtree(cfg["data_dir"] / "result" / "hic" / "bedpe", ignore_errors=True)
    (cfg["data_dir"] / "result" / "hic" / "bedpe").mkdir(parents=True, exist_ok=True)
    select_total_count = SelectTotalCount(cfg)
    df_merge = get_merge_bam(cfg)
    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        pairs2bedpe(
            pairs_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "pairs"
            / f"{exp}_{protein}_{treat}.pairs",
            bedpe_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_{treat}.bedpe",
            total_count=select_total_count(exp, protein, treat),
        )


def diff_bedpe(cfg: dict) -> None:
    df_merge = (
        get_merge_bam(cfg).query("treat.str.endswith('treat')").reset_index(drop=True)
    )
    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        diff = treat2diff(treat)
        treat_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_{treat}.bedpe"
        )
        control_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{'_'.join(map_to_wild_type_merge(exp, protein, treat))}.bedpe"
        )

        df = substract_bedpe(treat_file, control_file)

        for filter, direction in [("score > 0", "up"), ("score < 0", "down")]:
            df.query(filter).assign(score=lambda df: df["score"].abs()).to_csv(
                (
                    cfg["data_dir"]
                    / "result"
                    / "hic"
                    / "bedpe"
                    / f"{exp}_{protein}_{diff}.{direction}.bedpe"
                ),
                sep="\t",
                index=False,
                header=False,
            )


def sum_bedpe(cfg: dict) -> None:
    df_merge = get_merge_bam(cfg)
    select_total_count = SelectTotalCount(cfg)
    for exp, treat in (
        df_merge[["exp", "treat"]].drop_duplicates().itertuples(index=False)
    ):
        df_merge_slice = df_merge.query("exp == @exp and treat == @treat").reset_index(
            drop=True
        )

        bedpe_files = f"{exp}_" + df_merge_slice["protein"] + f"_{treat}.bedpe"
        total_counts = [
            select_total_count(exp, protein, treat)
            for protein in df_merge_slice["protein"]
        ]
        df_sum = summation_bedpe(
            [
                cfg["data_dir"] / "result" / "hic" / "bedpe" / bedpe_file
                for bedpe_file in bedpe_files
            ],
            total_counts,
        )

        df_sum.to_csv(
            cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp}_merge_{treat}.bedpe",
            sep="\t",
            index=False,
            header=False,
        )


def filter_non_cpcdh_junction(cfg: dict) -> None:
    (cfg["data_dir"] / "result" / "hic" / "bedpe" / "cpcdh").mkdir(
        exist_ok=True, parents=True
    )
    df_merge = get_merge_bam(cfg)
    for exp, treat in (
        df_merge[["exp", "treat"]].drop_duplicates().itertuples(index=False)
    ):
        assemble = treat2assemble(treat)
        df_intron = get_cpcdh_intron(
            cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
        )
        df_merge_slice = df_merge.query("exp == @exp and treat == @treat").reset_index(
            drop=True
        )
        for protein in ["merge"] + df_merge_slice["protein"].tolist():
            df = pd.read_csv(
                cfg["data_dir"]
                / "result"
                / "hic"
                / "bedpe"
                / f"{exp}_{protein}_{treat}.bedpe",
                sep="\t",
                names=[
                    "chrom1",
                    "start1",
                    "end1",
                    "chrom2",
                    "start2",
                    "end2",
                    "name",
                    "score",
                    "strand1",
                    "strand2",
                ],
            )

            df = (
                df
                .assign(**{
                    name: lambda df, start=start, end=end: (
                        (df["start1"] == start) & (df["start2"] == end)
                    )
                    for start, end, name in zip(
                        df_intron["start"],
                        df_intron["end"],
                        df_intron["name"],
                    )
                })
                .assign(
                    cpcdh=lambda df, df_intron=df_intron: df[
                        df_intron["name"].tolist()
                    ].any(axis=1)
                )
                .query("cpcdh")
                .reset_index(drop=True)[
                    [
                        "chrom1",
                        "start1",
                        "end1",
                        "chrom2",
                        "start2",
                        "end2",
                        "name",
                        "score",
                        "strand1",
                        "strand2",
                    ]
                ]
            )

            df.to_csv(
                cfg["data_dir"]
                / "result"
                / "hic"
                / "bedpe"
                / "cpcdh"
                / f"{exp}_{protein}_{treat}.bedpe",
                sep="\t",
                index=False,
                header=False,
            )

import os
import shutil

import pandas as pd
import py2bit

from .common import BlatSplice, get_cpcdh_intron, substract_bedpe, summation_bedpe
from .utils import (
    SelectTotalCount,
    clone2assemble,
    clone2treat,
    get_merge_bam,
    map_to_wild_type_merge,
    treat2assemble,
    treat2diff,
)


def get_pairs(cfg: dict) -> None:
    def splice_pair(ref_blocks: str) -> list:
        ref_blocks = ref_blocks.split(";")
        pairs = []
        for i in range(len(ref_blocks) - 1):
            block_chrom, block_start, block_end, block_strand = ref_blocks[i].split(":")
            block_start, block_end = int(block_start), int(block_end)
            next_block_chrom, next_block_start, next_block_end, next_block_strand = (
                ref_blocks[i + 1].split(":")
            )
            next_block_start, next_block_end = (
                int(next_block_start),
                int(next_block_end),
            )
            assert block_chrom == next_block_chrom, "inconsistent chrom"
            chrom1 = block_chrom
            chrom2 = next_block_chrom
            if block_start < next_block_start:
                pos1 = block_end + 1
                strand1 = block_strand
                pos2 = next_block_start + 1
                strand2 = next_block_strand
            else:
                pos1 = next_block_end + 1
                strand1 = next_block_strand
                pos2 = block_start + 1
                strand2 = block_strand

            pairs.append(f"{chrom1}:{pos1}:{strand1}:{chrom2}:{pos2}:{strand2}")

        return pairs

    df = (
        pd
        .read_feather(cfg["data_dir"] / "result" / "reads.feather")
        .query("ref_blocks.str.contains(';')")
        .reset_index(drop=True)
        .assign(
            treat=lambda df: df["clone"].map(clone2treat),
            exp_protein_treat=lambda df: (
                df["exp"] + "_" + df["protein"] + "_" + df["treat"]
            ),
            pair=lambda df: df["ref_blocks"].map(splice_pair),
        )[["exp_protein_treat", "query_name", "pair"]]
        .explode("pair", ignore_index=True)
    )

    df = pd.concat(
        [
            df[["exp_protein_treat", "query_name"]].rename(
                columns={"query_name": "readID"}
            ),
            df["pair"]
            .str.split(":", expand=True)
            .rename(
                columns={
                    0: "chrom1",
                    1: "pos1",
                    2: "strand1",
                    3: "chrom2",
                    4: "pos2",
                    5: "strand2",
                }
            )[["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]],
        ],
        axis=1,
    )

    (cfg["data_dir"] / "result" / "hic" / "pairs").mkdir(exist_ok=True, parents=True)
    for exp_protein_treat in df["exp_protein_treat"].unique():
        with open(
            cfg["data_dir"] / "result" / "hic" / "pairs" / f"{exp_protein_treat}.pairs",
            "w",
        ) as fd:
            fd.write("## pairs format v1.0\n")
            df.query("exp_protein_treat == @exp_protein_treat")[
                ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
            ].to_csv(fd, sep="\t", header=False, index=False)


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
                f"{chrom}:{chromStart}:{chromEnd}:{sourceChrom}:{sourceStart}:{sourceEnd}:{sourceStrand}:{targetChrom}:{targetStart}:{targetEnd}:{targetStrand}:{joint_block}"
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
                }
            ),
        ],
        axis=1,
    )

    blat_splice = BlatSplice(cfg)
    df_pidents = []
    for assemble in df["assemble"].unique():
        df_slice = df.query("assemble == @assemble")[
            ["uid", "joint_block"]
        ].reset_index(drop=True)
        blat_input = "\n".join(">" + df_slice["uid"] + "\n" + df_slice["joint_block"])
        blat_result = blat_splice(blat_input, assemble)
        df_pident = blat_result.groupby("qseqid", as_index=False)["pident"].max()
        df_pidents.append(df_pidents)
        breakpoint()

    df = df.assign(
        color=0,
        value=1.0,
        sourceName=".",
        targetName=".",
        joint_block_assemble=lambda df: df["joint_block"] + "|" + df["assemble"],
        pident=lambda df: df["joint_block_assemble"].map(
            lambda joint_block_assemble: blat_splice(*joint_block_assemble.split("|"))
        ),
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
        ]
    ].sort_values(by=["chrom", "chromStart"], ignore_index=True)

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


def pairs_to_bedpe(cfg: dict) -> None:
    shutil.rmtree(cfg["data_dir"] / "result" / "hic" / "bedpe", ignore_errors=True)
    (cfg["data_dir"] / "result" / "hic" / "bedpe").mkdir(parents=True, exist_ok=True)
    select_total_count = SelectTotalCount(cfg)
    for pairs_file in os.listdir(cfg["data_dir"] / "result" / "hic" / "pairs"):
        exp, protein, treat = pairs_file.removesuffix(".pairs").split("_")
        assemble = treat2assemble(treat)
        total_count = select_total_count(exp, protein, treat)

        pairs_file = cfg["data_dir"] / "result" / "hic" / "pairs" / pairs_file

        df = pd.read_csv(
            pairs_file,
            sep="\t",
            skiprows=1,
            names=[
                "readID",
                "chrom1",
                "pos1",
                "chrom2",
                "pos2",
                "strand1",
                "strand2",
            ],
        )

        df = (
            df
            .rename(
                columns={
                    "pos1": "end1",
                    "pos2": "end2",
                }
            )
            .assign(
                start1=lambda df: df["end1"] - 1,
                start2=lambda df: df["end2"] - 1,
            )
            .groupby(
                [
                    "chrom1",
                    "start1",
                    "end1",
                    "chrom2",
                    "start2",
                    "end2",
                ],
                as_index=False,
            )
            .agg(
                name=pd.NamedAgg("readID", lambda se: "|".join(se.tolist())),
                score=pd.NamedAgg("readID", "count"),
            )
            .assign(
                score=lambda df, total_count=total_count: (
                    df["score"] / total_count * 1_000_000
                )
            )[
                [
                    "chrom1",
                    "start1",
                    "end1",
                    "chrom2",
                    "start2",
                    "end2",
                    "name",
                    "score",
                ]
            ]
        )

        with py2bit.open(cfg[assemble]["2bit"]) as tb:
            df = df.assign(
                donor=lambda df: [
                    tb.sequence(chrom1, start1, start1 + 2)
                    for chrom1, start1 in zip(df["chrom1"], df["start1"])
                ],
                acceptor=lambda df: [
                    tb.sequence(chrom2, start2 - 2, start2)
                    for chrom2, start2 in zip(df["chrom2"], df["start2"])
                ],
            )

        df["strand1"] = "*"
        df["strand1"] = df["strand1"].where(
            (df["donor"] != "GT") | (df["acceptor"] != "AG"), "+"
        )
        df["strand1"] = df["strand1"].where(
            (df["donor"] != "CT") | (df["acceptor"] != "AC"), "-"
        )
        df["strand2"] = df["strand1"]

        df[
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
        ].to_csv(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / pairs_file.with_suffix(".bedpe").name,
            sep="\t",
            header=False,
            index=False,
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

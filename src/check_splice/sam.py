import os
import pathlib

import pandas as pd
import pysam
import sh

from . import check
from .common import (
    BlatSplice,
    ParseSamRead,
    filter_bam_reads,
    get_cpcdh_intron,
    is_mapped_primary_first,
)
from .utils import (
    clone2assemble,
    clone2treat,
    get_merge_bam,
    get_sample_bam,
    treat2assemble,
)


def parse_strand_sensitive_bam(cfg: dict):
    parse_sam_read = ParseSamRead()
    yield "exp,protein,clone,rep,query_name,query,is_forward,is_read1,is_qcfail,is_duplicate,mapping_quality,ref_block_chrom,ref_block_start,ref_block_end,ref_block_strand,query_block_start,query_block_end,align_string"
    for exp, protein, clone, rep, bamfile in get_sample_bam(cfg).itertuples(
        index=False
    ):
        bamfile = pathlib.Path(bamfile)
        assemble = clone2assemble(clone)
        for read in filter_bam_reads(
            bamfile,
            cfg[assemble]["chrom"],
            cfg[assemble]["start"],
            cfg[assemble]["end"],
        ):
            for (
                ref_block_chrom,
                ref_block_start,
                ref_block_end,
                ref_block_strand,
                query_block_start,
                query_block_end,
                align_string,
                query,
            ) in parse_sam_read.parse_read(read):
                if read.is_read1:
                    (
                        ref_block_chrom,
                        ref_block_start,
                        ref_block_end,
                        ref_block_strand,
                        query_block_start,
                        query_block_end,
                        align_string,
                        query,
                    ) = parse_sam_read.flip_read(
                        ref_block_chrom,
                        ref_block_start,
                        ref_block_end,
                        ref_block_strand,
                        query_block_start,
                        query_block_end,
                        align_string,
                        query,
                        read.query_length,
                    )

                yield f"{exp},{protein},{clone},{rep},{read.query_name},{query},{read.is_forward},{read.is_read1},{read.is_qcfail},{read.is_duplicate},{read.mapping_quality},{ref_block_chrom},{ref_block_start},{ref_block_end},{ref_block_strand},{query_block_start},{query_block_end},{align_string}"


def group_read_blocks(cfg: dict):
    df = pd.read_csv(cfg["data_dir"] / "result" / "reads.csv", header=0)
    df = (
        df
        .assign(
            ref_block=lambda df: (
                df["ref_block_chrom"]
                + ":"
                + df["ref_block_start"].astype(str)
                + ":"
                + df["ref_block_end"].astype(str)
                + ":"
                + df["ref_block_strand"]
            ),
            query_block=lambda df: (
                df["query_block_start"].astype(str)
                + ":"
                + df["query_block_end"].astype(str)
            ),
        )
        .sort_values(by="query_block_start")
        .groupby(
            by=[
                "exp",
                "protein",
                "clone",
                "rep",
                "query_name",
                "is_read1",
            ],
            as_index=False,
            sort=True,
        )
        .agg(
            query=pd.NamedAgg(column="query", aggfunc="first"),
            is_forward=pd.NamedAgg(column="is_forward", aggfunc="first"),
            is_qcfail=pd.NamedAgg(column="is_qcfail", aggfunc="first"),
            is_duplicate=pd.NamedAgg(column="is_duplicate", aggfunc="first"),
            mapping_quality=pd.NamedAgg(column="mapping_quality", aggfunc="first"),
            ref_blocks=pd.NamedAgg(column="ref_block", aggfunc=";".join),
            query_blocks=pd.NamedAgg(column="query_block", aggfunc=";".join),
            align_strings=pd.NamedAgg(column="align_string", aggfunc=";".join),
        )
    )
    df.to_feather(cfg["data_dir"] / "result" / "reads.feather")


class ExpandSplice:
    def __init__(self) -> None:
        pass

    def parse_blocks(self, ref_query_blocks: str) -> list:
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

            interact.append(
                f"{ref_block_chrom}:{ref_block_start}:{ref_block_end}:{ref_block_strand}:{next_ref_block_chrom}:{next_ref_block_start}:{next_ref_block_end}:{next_ref_block_strand}:{query[query_block_start:next_query_block_end]}:{i}"
            )

        return interact

    def __call__(self, cfg: dict):
        df = (
            pd
            .read_feather(cfg["data_dir"] / "result" / "reads.feather")
            .query("ref_blocks.str.contains(';')")
            .reset_index(drop=True)
            .assign(
                ref_query_blocks=lambda df: (
                    df["ref_blocks"] + "|" + df["query_blocks"] + "|" + df["query"]
                ),
                parse=lambda df: df["ref_query_blocks"].map(self.parse_blocks),
            )[["exp", "protein", "clone", "rep", "query_name", "is_read1", "parse"]]
            .explode("parse", ignore_index=True)
        )

        df = pd.concat(
            [
                df,
                df["parse"]
                .str.split(":", expand=True)
                .rename(
                    columns={
                        0: "chrom1",
                        1: "start1",
                        2: "end1",
                        3: "strand1",
                        4: "chrom2",
                        5: "start2",
                        6: "end2",
                        7: "strand2",
                        8: "joint_block",
                        9: "bidx",
                    }
                )
                .astype({
                    "start1": int,
                    "end1": int,
                    "start2": int,
                    "end2": int,
                    "bidx": int,
                }),
            ],
            axis=1,
        )

        df = df.assign(
            assemble=lambda df: df["clone"].map(clone2assemble),
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
                + "_"
                + df["bidx"].astype(str)
            ),
        )

        blat_splice = BlatSplice(cfg)
        blat_results = []
        for assemble in df["assemble"].unique():
            df_slice = df.query("assemble == @assemble").reset_index(drop=True)
            blat_input = "\n".join(
                ">" + df_slice["uid"] + "\n" + df_slice["joint_block"]
            )
            blat_result = blat_splice(blat_input, assemble)
            blat_results.append(blat_result)

        df = df.merge(
            pd.concat(blat_results, ignore_index=True),
            how="left",
            left_on="uid",
            right_on="qseqid",
            validate="one_to_one",
        ).assign(
            match_base=lambda df: df["match_base"].fillna(0.0),
            detail=lambda df: df["detail"].fillna(""),
        )

        df[
            [
                "exp",
                "protein",
                "clone",
                "rep",
                "query_name",
                "is_read1",
                "bidx",
                "chrom1",
                "start1",
                "end1",
                "strand1",
                "chrom2",
                "start2",
                "end2",
                "strand2",
                "joint_block",
                "match_base",
                "detail",
            ]
        ].to_csv(cfg["data_dir"] / "result" / "expand_splice.csv", index=False)


def merge_bam(cfg: dict) -> None:
    (cfg["data_dir"] / "bam" / "merge").mkdir(exist_ok=True, parents=True)
    samtools = sh.Command("samtools")
    df = get_sample_bam(cfg).assign(
        treat=lambda df: df["clone"].map(clone2treat),
        assemble=lambda df: df["clone"].map(clone2assemble),
    )
    for assemble in ["hg19", "mm10"]:
        for exp in ["total", "rna", "pro", "clip"]:
            for protein in ["WT", "NP220", "MPP8", "PPHLN1", "TASOR"]:
                for treat in ["control", "treat"]:
                    if assemble == "mm10":
                        treat = f"mm{treat}"

                    bamfiles = df.query(
                        "assemble == @assemble and exp == @exp and protein == @protein and treat == @treat"
                    )["file"].tolist()
                    if not bamfiles:
                        continue

                    merge_bam = (
                        cfg["data_dir"]
                        / "bam"
                        / "merge"
                        / f"{exp}_{protein}_{treat}.bam"
                    )
                    print(merge_bam, bamfiles)

                    samtools(
                        "merge",
                        "-f",
                        "-o",
                        os.fspath(merge_bam),
                        *bamfiles,
                    )
                    samtools("index", os.fspath(merge_bam))


def filter_bam(cfg: dict) -> None:
    (cfg["data_dir"] / "bam" / "merge" / "precursor").mkdir(parents=True, exist_ok=True)
    (cfg["data_dir"] / "bam" / "merge" / "splice").mkdir(parents=True, exist_ok=True)
    samtools = sh.Command("samtools")
    parse_sam_read = ParseSamRead()
    for exp, protein, treat, bamfile in get_merge_bam(cfg).itertuples(index=False):
        bamfile = pathlib.Path(bamfile)
        assemble = treat2assemble(treat)
        df_intron = get_cpcdh_intron(
            cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
        )
        intervals = df_intron.melt(
            id_vars=["chrom", "name"],
            value_vars=["start", "end"],
            var_name="se",
            value_name="pos",
        ).assign(
            start=lambda df: df["pos"] - cfg["cover_threshold"],
            end=lambda df: df["pos"] + cfg["cover_threshold"],
        )["chrom", "start", "end"]

        with pysam.AlignmentFile(bamfile, "rb") as infile:
            precursor_forward_bam_file = (
                bamfile.with_name("precursor") / bamfile.with_suffix(".f.bam").name
            )
            precursor_reverse_bam_file = (
                bamfile.with_name("precursor") / bamfile.with_suffix(".r.bam").name
            )
            splice_bam_file = bamfile.with_name("splice") / bamfile.name
            with (
                pysam.AlignmentFile(
                    precursor_forward_bam_file, "wb", template=infile
                ) as pfb,
                pysam.AlignmentFile(
                    precursor_reverse_bam_file, "wb", template=infile
                ) as prb,
                pysam.AlignmentFile(splice_bam_file, "wb", template=infile) as sb,
            ):
                for read in infile.fetch(
                    cfg[assemble]["chrom"],
                    cfg[assemble]["start"],
                    cfg[assemble]["end"],
                ):
                    if not is_mapped_primary_first(read):
                        continue

                    if check.any_cover_any_nostrand(
                        [
                            (
                                ref_block_chrom,
                                ref_block_start,
                                ref_block_end,
                                ref_block_strand,
                            )
                            for ref_block_chrom, ref_block_start, ref_block_end, ref_block_strand, query_block_start, query_block_end, align_string in parse_sam_read.parse_read()
                        ],
                        intervals,
                    ):
                        if (
                            read.is_read1
                            and read.is_reverse
                            or read.is_read2
                            and read.is_forward
                        ):
                            pfb.write(read)
                        else:
                            prb.write(read)

                    if "N" in read.cigarstring or read.has_tag("SA"):
                        sb.write(read)

        samtools("index", os.fspath(precursor_forward_bam_file))
        samtools("index", os.fspath(precursor_reverse_bam_file))
        samtools("index", os.fspath(splice_bam_file))

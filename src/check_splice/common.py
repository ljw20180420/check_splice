import os
import re
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pysam
import sh


def get_bam_read_count(bamfile: os.PathLike) -> int:
    with pysam.AlignmentFile(os.fspath(bamfile)) as bam:
        total_count = sum(
            1
            for read in bam
            if not read.is_secondary and read.is_mapped and not read.is_supplementary
        )

    return total_count


def prepare_gene_bed12(
    gtffile: os.PathLike, outfile: os.PathLike, addtional_filter: Callable
) -> None:
    gtffile = Path(os.fspath(gtffile))
    outfile = Path(os.fspath(outfile))

    gtfToGenePred = sh.Command("gtfToGenePred")
    genePredToBigGenePred = sh.Command("genePredToBigGenePred")
    gtfToGenePred(
        "-genePredExt",
        os.fspath(gtffile),
        os.fspath(gtffile.with_suffix(".gp")),
    )
    df_gp = pd.read_csv(gtffile.with_suffix(".gp"), sep="\t", header=None)
    df_gp[[11] + list(range(1, 15))].to_csv(
        gtffile.with_suffix(".gp"),
        sep="\t",
        header=False,
        index=False,
    )
    genePredToBigGenePred(
        os.fspath(gtffile.with_suffix(".gp")),
        os.fspath(gtffile.with_suffix(".bgp")),
    )
    df_bgp = (
        pd
        .read_csv(gtffile.with_suffix(".bgp"), sep="\t", header=None)[list(range(12))]
        .rename(
            columns={
                0: "chrom",
                1: "chromStart",
                2: "chromEnd",
                3: "name",
                4: "score",
                5: "strand",
                6: "thickStart",
                7: "thickEnd",
                8: "itemRgb",
                9: "blockCount",
                10: "blockSizes",
                11: "blockStarts",
            }
        )
        .sort_values(by=["chrom", "chromStart"], ignore_index=True)
    )

    df_bgp = (
        df_bgp
        .query("not name.str.lower().str.startswith('pcdha') or blockCount == 4")
        .query(
            "not name.str.lower().str.startswith('pcdhb') or blockCount == 1 or name.str.lower() == 'pcdhb9'"
        )
        .query("not name.str.lower().str.startswith('pcdhb@')")
        .query("not name.str.lower().str.startswith('pcdhg') or blockCount == 4")
        .reset_index(drop=True)
    )

    df_bgp = addtional_filter(df_bgp)

    df_bgp.to_csv(
        outfile,
        sep="\t",
        header=False,
        index=False,
    )
    outfile.with_suffix(".bed.bgz").unlink(missing_ok=True)
    outfile.with_suffix(".bed.bgz.tbi").unlink(missing_ok=True)


def get_cpcdh_exon(gtffile: os.PathLike, chrom: str) -> pd.DataFrame:
    df = pd.read_csv(
        gtffile,
        sep="\t",
        names=[
            "chrom",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ],
    )
    df = df.query(
        "chrom == @chrom and attributes.str.lower().str.contains(r'pcdh[abg][abc]?[0-9]{1,2}') and (feature=='exon' or feature=='CDS')"
    ).reset_index(drop=True)
    attributes = df["attributes"].str.split(expand=True)
    df = df.assign(
        start=lambda df: df["start"] - 1,
        name=attributes[9].str.strip('";'),
        transcript_id=attributes[3].str.strip('";'),
        exon_number=attributes[5].str.strip('";').astype(int),
        total_exon_number=lambda df: df.groupby("transcript_id")[
            "exon_number"
        ].transform(max),
    ).drop(
        columns=[
            "source",
            "score",
            "frame",
            "attributes",
        ]
    )
    df = (
        df
        .query(
            "(exon_number == 1 and (total_exon_number == 4 or name.str.lower().str.startswith('pcdhb'))) or (exon_number > 1 and (name.str.lower() == 'pcdha1' or name.str.lower() == 'pcdhga1'))"
        )
        .reset_index(drop=True)
        .assign(
            repeat=lambda df: df.groupby(["feature", "name", "exon_number"])[
                "name"
            ].transform("count"),
        )
        .query("repeat == 1 or transcript_id.str.contains(r'(?:^NM_018|NM_002)')")
        .reset_index(drop=True)
        .assign(
            name=lambda df: df.apply(
                lambda row: (
                    row["name"]
                    if row["exon_number"] == 1
                    else f"ace{row['exon_number'] - 1}"
                    if row["name"].lower() == "pcdha1"
                    else f"gce{row['exon_number'] - 1}"
                ),
                axis=1,
            )
        )
    )
    df = df.pivot_table(
        values=["start", "end"],
        index=["chrom", "strand", "name", "transcript_id"],
        columns="feature",
    )
    df.columns = df.columns.to_flat_index().map(lambda tp: f"{tp[1]}_{tp[0]}")
    df = (
        df
        .reset_index()
        .rename(columns={"exon_start": "start", "exon_end": "end"})
        .assign(
            CDS_start=lambda df: df["CDS_start"].fillna("."),
            CDS_end=lambda df: df["CDS_end"].fillna("."),
        )
        .assign(score=".")[
            [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "CDS_start",
                "CDS_end",
                "transcript_id",
            ]
        ]
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["start", "end"], ignore_index=True)
    )

    return df


def get_cpcdh_intron(cpcdh_csv: os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(cpcdh_csv, header=0)

    intron_names = []
    intron_starts = []
    intron_ends = []

    intron_end = df.query("name == 'ace1'")["start"].item()
    for name, intron_start in df.query("name.str.lower().str.startswith('pcdha')")[
        ["name", "end"]
    ].itertuples(index=False):
        intron_names.append(f"{name}_ace1")
        intron_starts.append(intron_start)
        intron_ends.append(intron_end)

    intron_names.append("ace1_ace2")
    intron_starts.append(df.query("name == 'ace1'")["end"].item())
    intron_ends.append(df.query("name == 'ace2'")["start"].item())

    intron_names.append("ace2_ace3")
    intron_starts.append(df.query("name == 'ace2'")["end"].item())
    intron_ends.append(df.query("name == 'ace3'")["start"].item())

    intron_end = df.query("name == 'gce1'")["start"].item()
    for name, intron_start in df.query("name.str.lower().str.startswith('pcdhg')")[
        ["name", "end"]
    ].itertuples(index=False):
        intron_names.append(f"{name}_gce1")
        intron_starts.append(intron_start)
        intron_ends.append(intron_end)

    intron_names.append("gce1_gce2")
    intron_starts.append(df.query("name == 'gce1'")["end"].item())
    intron_ends.append(df.query("name == 'gce2'")["start"].item())

    intron_names.append("gce2_gce3")
    intron_starts.append(df.query("name == 'gce2'")["end"].item())
    intron_ends.append(df.query("name == 'gce3'")["start"].item())

    return (
        pd
        .DataFrame({
            "start": intron_starts,
            "end": intron_ends,
            "name": intron_names,
        })
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["start"], ignore_index=True)
    )


def filter_bam_reads(bamfile: os.PathLike, chrom: str, start: int, end: int):
    with pysam.AlignmentFile(os.fspath(bamfile)) as fd:
        for read in fd.fetch(
            contig=chrom,
            start=start,
            end=end,
        ):
            if read.is_secondary:
                continue
            if not read.is_mapped:
                continue
            if read.is_supplementary:
                continue

            yield read


class ParseSamRead:
    def __init__(self):
        self.cigar_parser = re.compile(r"(\d+)([MIDNSHP=XB])")
        self.align_sting_parser = re.compile(r"\d+[MID]")

    def flip_to_query_raw_strand(
        self,
        query_block_start: int,
        query_block_end: int,
        align_string: list[str],
        strand: str,
        query_length: int,
    ):
        if strand == "+":
            return (
                query_block_start,
                query_block_end,
                "".join(align_string),
            )
        else:
            # strand == "-"
            return (
                query_length - query_block_end,
                query_length - query_block_start,
                "".join(reversed(align_string)),
            )

    def parse_cigar(self, start: int, cigarstring: str, strand: str, query_length: int):
        ref_current_pos = start
        ref_block_start = start
        query_current_pos = 0
        length_ops = [
            (int(length), op) for length, op in self.cigar_parser.findall(cigarstring)
        ]
        if length_ops[0][1] == "S":
            query_current_pos += length_ops[0][0]
            length_ops = length_ops[1:]
        query_block_start = query_current_pos

        align_string = []
        for length, op in length_ops:
            if op in ("=", "X"):
                op = "M"

            if op in ("M", "I", "D"):
                align_string.append(f"{length}{op}")
                if op in ("M", "D"):
                    ref_current_pos += length
                if op in ("M", "I"):
                    query_current_pos += length

            elif op == "N":
                yield (
                    ref_block_start,
                    ref_current_pos,
                    *self.flip_to_query_raw_strand(
                        query_block_start,
                        query_current_pos,
                        align_string,
                        strand,
                        query_length,
                    ),
                )

                ref_current_pos += length
                ref_block_start = ref_current_pos
                query_block_start = query_current_pos
                align_string = []

        yield (
            ref_block_start,
            ref_current_pos,
            *self.flip_to_query_raw_strand(
                query_block_start,
                query_current_pos,
                align_string,
                strand,
                query_length,
            ),
        )

    def parse_sa(self, sa_tag: str):
        for alignment_str in sa_tag.split(";"):
            if not alignment_str:
                continue

            chrom, start, strand, cigar, _ = alignment_str.split(",")
            start = int(start) - 1

            yield chrom, start, strand, cigar

    def parse_read(self, read: pysam.AlignedSegment):
        chroms = [read.reference_name]
        starts = [read.reference_start]
        strands = ["+" if read.is_forward else "-"]
        cigars = [read.cigarstring]
        if read.has_tag("SA"):
            for chrom, start, strand, cigar in self.parse_sa(read.get_tag("SA")):
                chroms.append(chrom)
                starts.append(start)
                strands.append(strand)
                cigars.append(cigar)

        for chrom, start, strand, cigar in zip(chroms, starts, strands, cigars):
            for (
                ref_block_start,
                ref_block_end,
                query_block_start,
                query_block_end,
                align_string,
            ) in self.parse_cigar(start, cigar, strand, read.query_length):
                yield (
                    chrom,
                    ref_block_start,
                    ref_block_end,
                    strand,
                    query_block_start,
                    query_block_end,
                    align_string,
                )

    def flip_read(
        self,
        ref_block_chrom: str,
        ref_block_start: int,
        ref_block_end: int,
        ref_block_strand: str,
        query_block_start: int,
        query_block_end: int,
        align_string: str,
        query_length: int,
    ):
        return (
            ref_block_chrom,
            ref_block_start,
            ref_block_end,
            "+" if ref_block_strand == "-" else "-",
            query_length - query_block_end,
            query_length - query_block_start,
            "".join(reversed(self.align_sting_parser.findall(align_string))),
        )

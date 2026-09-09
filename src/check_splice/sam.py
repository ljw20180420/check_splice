import os
import re

import pandas as pd
import pyBigWig
import pysam
import sh

from .utils import get_precursor_pos


def parse_cigar(start: int, cigarstring: str):
    blocks = []
    current_pos = start
    block_start = start
    pattern = re.compile(r"(\d+)([MIDNSHP=XB])")
    for length, op in pattern.findall(cigarstring):
        length = int(length)
        if op in ("M", "D", "=", "X"):  # Operators that consume reference genome space
            current_pos += length
        elif op == "N":  # Intron / Reference Skip (N)
            # End the current block before the intron starts
            if current_pos > block_start:
                blocks.append((block_start, current_pos))
            # Skip past the intron region
            current_pos += length
            # Set the start of the next block to the end of the intron
            block_start = current_pos
        elif op not in ("I", "S", "H"):
            # Insertions and clips do not move the reference cursor
            raise ValueError("unknown cigar operation")

    # Append the final block after parsing the last CIGAR operation
    if current_pos > block_start:
        blocks.append((block_start, current_pos))

    return blocks


def parse_sa(sa_tag: str):
    for alignment_str in sa_tag.split(";"):
        if not alignment_str:
            continue

        chrom, start, strand, cigar, _ = alignment_str.split(",")
        start = int(start) - 1

        yield chrom, start, strand, cigar


def parse_block_without_flip(read: pysam.AlignedSegment):
    chroms = [read.reference_name]
    starts = [read.reference_start]
    strands = ["+" if read.is_forward else "-"]
    cigars = [read.cigarstring]
    if read.has_tag("SA"):
        for chrom, start, strand, cigar in parse_sa(read.get_tag("SA")):
            chroms.append(chrom)
            starts.append(start)
            strands.append(strand)
            cigars.append(cigar)

    for chrom, start, strand, cigar in zip(chroms, starts, strands, cigars):
        blocks = parse_cigar(start, cigar)
        if strand == "+":
            for block_start, block_end in blocks:
                yield chrom, block_start, block_end, strand
        else:
            for block_start, block_end in reversed(blocks):
                yield chrom, block_start, block_end, strand


def parse_block_with_flip(read: pysam.AlignedSegment, flip: str) -> list:
    blocks = list(parse_block_without_flip(read))
    assert flip in ("R1", "R2"), "flip must be either 'R1' or 'R2'"
    if flip == "R2" and read.is_read1 or flip == "R1" and read.is_read2:
        return blocks

    flip_blocks = []
    for chrom, block_start, block_end, strand in reversed(blocks):
        flip_blocks.append((
            chrom,
            block_start,
            block_end,
            "+" if strand == "-" else "-",
        ))

    return flip_blocks


def filter_reads(samfile: os.PathLike, chrom: str, start: int, end: int):
    with pysam.AlignmentFile(os.fspath(samfile)) as fd:
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


def merge_bam(cfg: dict) -> None:
    (cfg["data_dir"] / "bam" / "merge").mkdir(exist_ok=True, parents=True)
    samtools = sh.Command("samtools")
    for exp in ["total", "rna", "pro", "clip"]:
        for protein in ["WT", "NP220", "MPP8", "PPHLN1", "TASOR"]:
            for wt in [True, False]:
                if wt:
                    treat = "control"
                else:
                    if exp != "clip":
                        treat = "delta"
                    else:
                        treat = "tag"

                bam_files = []
                for bam_file in os.listdir(cfg["data_dir"] / "bam"):
                    if not bam_file.endswith(".bam"):
                        continue
                    exp_, protein_, clone_, _ = bam_file.split("_", 3)
                    if exp_ != exp or protein_ != protein:
                        continue
                    if wt != clone_.startswith("WT"):
                        continue

                    bam_file = cfg["data_dir"] / "bam" / bam_file
                    bam_files.append(os.fspath(bam_file))

                if not bam_files:
                    continue

                merge_bam = (
                    cfg["data_dir"] / "bam" / "merge" / f"{exp}_{protein}_{treat}.bam"
                )
                samtools(
                    "merge",
                    "-f",
                    "-o",
                    os.fspath(merge_bam),
                    *bam_files,
                )
                samtools("index", os.fspath(merge_bam))


def filter_precursor_bam(cfg: dict, bam_file: os.PathLike, strand: str) -> None:
    df_se = get_precursor_pos(cfg)
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        filtered_bam_file = (
            bam_file.with_name("precursor")
            / bam_file.with_suffix(f".{strand}.bam").name
        )
        with pysam.AlignmentFile(filtered_bam_file, "wb", template=infile) as outfile:
            for read in infile.fetch(cfg["chrom"], cfg["start"], cfg["end"]):
                if read.is_secondary:
                    continue
                if not read.is_mapped:
                    continue
                if read.is_supplementary:
                    continue

                read_strand = (
                    "f"
                    if read.is_read1
                    and read.is_reverse
                    or read.is_read2
                    and read.is_forward
                    else "r"
                )
                if read_strand != strand:
                    continue

                for (
                    block_chrom,
                    block_start,
                    block_end,
                    block_strand,
                ) in parse_block_without_flip(read):
                    cover_up = (df_se["pos"] + cfg["cover_threshold"]).between(
                        block_start, block_end
                    )
                    cover_down = (df_se["pos"] - cfg["cover_threshold"]).between(
                        block_start, block_end
                    )

                    if (cover_up & cover_down).any():
                        outfile.write(read)
                        break

    samtools = sh.Command("samtools")
    samtools(
        "index",
        os.fspath(filtered_bam_file),
    )


def filter_splice_bam(cfg: dict, bam_file: os.PathLike) -> None:
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        filtered_bam_file = bam_file.with_name("splice") / bam_file.name
        with pysam.AlignmentFile(filtered_bam_file, "wb", template=infile) as outfile:
            for read in infile.fetch(cfg["chrom"], cfg["start"], cfg["end"]):
                if read.is_secondary:
                    continue
                if not read.is_mapped:
                    continue
                if read.is_supplementary:
                    continue

                for i, (
                    block_chrom,
                    block_start,
                    block_end,
                    block_strand,
                ) in enumerate(parse_block_without_flip(read)):
                    if i > 0:
                        outfile.write(read)
                        break

    samtools = sh.Command("samtools")
    samtools(
        "index",
        os.fspath(filtered_bam_file),
    )


def filter_bam_all(cfg: dict) -> None:
    (cfg["data_dir"] / "bam" / "merge" / "precursor").mkdir(parents=True, exist_ok=True)
    (cfg["data_dir"] / "bam" / "merge" / "splice").mkdir(parents=True, exist_ok=True)
    for bam_file in os.listdir(cfg["data_dir"] / "bam" / "merge"):
        if not bam_file.endswith(".bam"):
            continue

        bam_file = cfg["data_dir"] / "bam" / "merge" / bam_file

        for strand in ["f", "r"]:
            filter_precursor_bam(cfg, bam_file, strand)

        filter_splice_bam(cfg, bam_file)

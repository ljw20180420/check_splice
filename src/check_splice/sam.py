import os
import re

import pandas as pd
import pysam
import sh

from .utils import get_precursor_pos


def parse_cigar(start: int, cigarstring: str, strand: str, query_length: int):
    blocks = []
    query_blocks = []
    align_strings = []
    current_pos = start
    block_start = start
    query_current_pos = 0
    pattern = re.compile(r"(\d+)([MIDNSHP=XB])")
    length_ops = [(int(length), op) for length, op in pattern.findall(cigarstring)]
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
                current_pos += length
            if op in ("M", "I"):
                query_current_pos += length

        elif op == "N":  # Intron / Reference Skip (N)
            # End the current block before the intron starts
            blocks.append((block_start, current_pos))
            query_blocks.append((query_block_start, query_current_pos))
            align_strings.append(align_string)
            # Skip past the intron region
            current_pos += length
            # Set the start of the next block to the end of the intron
            block_start = current_pos
            query_block_start = query_current_pos
            align_string = []

    blocks.append((block_start, current_pos))
    query_blocks.append((query_block_start, query_current_pos))
    align_strings.append("".join(align_string))

    if strand == "-":
        query_blocks = [
            (query_length - query_block_end, query_length - query_block_start)
            for query_block_start, query_block_end in query_blocks
        ]
        align_strings = [
            "".join(reversed(align_string)) for align_string in align_strings
        ]
    else:
        align_strings = ["".join(align_string) for align_string in align_strings]

    return blocks, query_blocks, align_strings


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

    all_chroms = []
    all_block_starts = []
    all_block_ends = []
    all_query_block_starts = []
    all_query_block_ends = []
    all_align_strings = []
    all_strands = []
    for chrom, start, strand, cigar in zip(chroms, starts, strands, cigars):
        blocks, query_blocks, align_strings = parse_cigar(
            start, cigar, strand, read.query_length
        )
        for (block_start, block_end), (
            query_block_start,
            query_block_end,
        ), align_string in zip(blocks, query_blocks, align_strings):
            all_chroms.append(chrom)
            all_block_starts.append(block_start)
            all_block_ends.append(block_end)
            all_query_block_starts.append(query_block_start)
            all_query_block_ends.append(query_block_end)
            all_align_strings.append(align_string)
            all_strands.append(strand)

    return pd.DataFrame({
        "chrom": all_chroms,
        "block_start": all_block_starts,
        "block_end": all_block_ends,
        "query_block_start": all_query_block_starts,
        "query_block_end": all_query_block_ends,
        "align_string": all_align_strings,
        "strand": all_strands,
    }).sort_values(by=["query_block_start"], ignore_index=True)


def parse_block_with_flip(
    read: pysam.AlignedSegment, flip: str
) -> tuple[list, list, list]:
    df = parse_block_without_flip(read)
    assert flip in ("R1", "R2"), "flip must be either 'R1' or 'R2'"
    if flip == "R1" and read.is_read1 or flip == "R2" and read.is_read2:
        pattern = re.compile(r"\d+[MID]")
        df = (
            df
            .assign(
                strand=lambda df: df["strand"].map({"+": "-", "-": "+"}),
            )
            .rename(
                columns={
                    "query_block_start": "query_block_start_old",
                    "query_block_end": "query_block_end_old",
                }
            )
            .assign(
                query_block_start=lambda df: (
                    read.query_length - df["query_block_end_old"]
                ),
                query_block_end=lambda df: (
                    read.query_length - df["query_block_start_old"]
                ),
            )
            .drop(columns=["query_block_start_old", "query_block_end_old"])
            .assign(
                align_string=lambda df: df["align_string"].map(
                    lambda align_string: "".join(
                        reversed(pattern.findall(align_string))
                    )
                )
            )
        )
        df = df[::-1].reset_index(drop=True)

    ref_blocks = list(
        zip(df["chrom"], df["block_start"], df["block_end"], df["strand"])
    )
    query_blocks = list(zip(df["query_block_start"], df["query_block_end"]))
    align_strings = list(df["align_string"])

    return ref_blocks, query_blocks, align_strings


def filter_reads(bamfile: os.PathLike, chrom: str, start: int, end: int):
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


def merge_bam(cfg: dict, assemble: str) -> None:
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
                if assemble == "mm10":
                    treat = f"mm{treat}"

                bam_files = []
                for bam_file in os.listdir(cfg["data_dir"] / "bam"):
                    if not bam_file.endswith(".bam"):
                        continue
                    exp_, protein_, clone_, _ = bam_file.split("_", 3)
                    if exp_ != exp or protein_ != protein:
                        continue
                    if (assemble == "mm10") != clone_.startswith("mm"):
                        continue
                    if (assemble != "mm10") and wt != clone_.startswith("WT"):
                        continue
                    if (assemble == "mm10") and wt != clone_.startswith("mmWT"):
                        continue

                    bam_file = cfg["data_dir"] / "bam" / bam_file
                    bam_files.append(os.fspath(bam_file))

                if not bam_files:
                    continue

                merge_bam = (
                    cfg["data_dir"] / "bam" / "merge" / f"{exp}_{protein}_{treat}.bam"
                )
                print(merge_bam, bam_files)

                samtools(
                    "merge",
                    "-f",
                    "-o",
                    os.fspath(merge_bam),
                    *bam_files,
                )
                samtools("index", os.fspath(merge_bam))


def filter_precursor_bam(cfg: dict, bam_file: os.PathLike, strand: str) -> None:
    exp, protein, treat = bam_file.name.removesuffix(".bam").split("_")
    assemble = "hg19" if not treat.startswith("mm") else "mm10"
    df_se = get_precursor_pos(cfg, assemble)
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        filtered_bam_file = (
            bam_file.with_name("precursor")
            / bam_file.with_suffix(f".{strand}.bam").name
        )
        with pysam.AlignmentFile(filtered_bam_file, "wb", template=infile) as outfile:
            for read in infile.fetch(
                cfg[assemble]["chrom"],
                cfg[assemble]["start"],
                cfg[assemble]["end"],
            ):
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
                    query_block_start,
                    query_block_end,
                    align_string,
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
    exp, protein, treat = bam_file.name.removesuffix(".bam").split("_")
    assemble = "hg19" if not treat.startswith("mm") else "mm10"
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        filtered_bam_file = bam_file.with_name("splice") / bam_file.name
        with pysam.AlignmentFile(filtered_bam_file, "wb", template=infile) as outfile:
            for read in infile.fetch(
                cfg[assemble]["chrom"], cfg[assemble]["start"], cfg[assemble]["end"]
            ):
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
                    query_block_start,
                    query_block_end,
                    align_string,
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

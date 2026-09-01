import os
import re

import pysam


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


def parse_block_without_flip_R2(read: pysam.AlignedSegment):
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
    blocks = list(parse_block_without_flip_R2(read))
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

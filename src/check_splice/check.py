import os

import pandas as pd


def map_start(blocks: list[tuple[str, int, int, str]]) -> int:
    _, block_start, block_end, block_strand = blocks[0]
    if block_strand == "+":
        return block_start
    else:
        return block_end


class Interval:
    def __init__(
        self, chrom: str, start: int, end: int, strand: str, name: str
    ) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name

    def connect(self, blocks: list[tuple[str, int, int, str]]) -> bool:
        for i in range(len(blocks) - 1):
            block_chrom, block_start, block_end, block_strand = blocks[i]
            next_block_chrom, next_block_start, next_block_end, next_block_strand = (
                blocks[i + 1]
            )
            if (
                self.chrom == block_chrom
                and self.chrom == next_block_chrom
                and self.strand == block_strand
                and self.strand == next_block_strand
            ):
                if self.strand == "+":
                    if self.start == block_end and self.end == next_block_start:
                        return True
                else:
                    if self.start == next_block_end and self.end == block_start:
                        return True
        return False

    def cover(self, blocks: list[tuple[str, int, int, str]]) -> bool:
        for block_chrom, block_start, block_end, block_strand in blocks:
            if (
                self.chrom == block_chrom
                and self.strand == block_strand
                and self.start >= block_start
                and self.end <= block_end
            ):
                return True
        return False

    def overlap(self, blocks: list[tuple[str, int, int, str]]) -> bool:
        for block_chrom, block_start, block_end, block_strand in blocks:
            if (
                self.chrom == block_chrom
                and self.strand == block_strand
                and min(self.end, block_end) > max(self.start, block_start)
            ):
                return True
        return False


class Intervals:
    def __init__(self, intervals: list[Interval]) -> None:
        self.intervals = intervals

    def __call__(
        self, info: dict, blocks: list[tuple[str, int, int, str]], method: str
    ) -> dict:
        for interval in self.intervals:
            func = getattr(interval, method)
            info[f"{method}_{interval.name}"] = func(blocks)

        return info


def all_intervals(cpcdh_file: os.PathLike, cover_threshold: int):
    df = pd.read_csv(cpcdh_file, header=0)
    introns = Intervals([
        Interval(chrom, start, end, strand, name)
        for chrom, start, end, strand, name in df.query(type="intron")[
            ["chrom", "start", "end", "strand", "name"]
        ].itertuples(index=False)
    ])
    intron_starts = Intervals([
        Interval(
            chrom,
            start - cover_threshold,
            start + cover_threshold,
            strand,
            f"{name}_start",
        )
        for chrom, start, end, strand, name in df.query(type="intron")[
            ["chrom", "start", "end", "strand", "name"]
        ].itertuples(index=False)
    ])
    intron_ends = Intervals([
        Interval(
            chrom,
            end - cover_threshold,
            end + cover_threshold,
            strand,
            f"{name}_end",
        )
        for chrom, start, end, strand, name in df.query(type="intron")[
            ["chrom", "start", "end", "strand", "name"]
        ].itertuples(index=False)
    ])

    return introns, intron_starts, intron_ends

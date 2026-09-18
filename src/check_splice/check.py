import pandas as pd


def cover(
    blocks: list[tuple[str, int, int, str]],
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> bool:
    for block_chrom, block_start, block_end, block_strand in blocks:
        if (
            chrom == block_chrom
            and strand == block_strand
            and start >= block_start
            and end <= block_end
        ):
            return True
    return False


def connect(
    blocks: list[tuple[str, int, int, str]],
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> bool:
    for i in range(len(blocks) - 1):
        block_chrom, block_start, block_end, block_strand = blocks[i]
        next_block_chrom, next_block_start, next_block_end, next_block_strand = blocks[
            i + 1
        ]
        if (
            chrom == block_chrom
            and chrom == next_block_chrom
            and strand == block_strand
            and strand == next_block_strand
        ):
            if strand == "+":
                if start == block_end and end == next_block_start:
                    return True
            else:
                if start == next_block_end and end == block_start:
                    return True
    return False


def start(blocks: list[tuple[str, int, int, str]]) -> int:
    chrom, block_start, block_end, strand = blocks[0]
    if strand == "+":
        return block_start
    else:
        # strand = "-"
        return block_end


def any_cover_any_nostrand(
    blocks: list[tuple[str, int, int, str]], intervals: pd.DataFrame
) -> bool:
    for chrom, start, end, strand in blocks:
        intervals_splice = intervals.query("chrom == @chrom").reset_index(drop=True)
        cover_up = (intervals_splice["start"]).between(start, end)
        cover_down = (intervals_splice["end"]).between(start, end)

        return (cover_up & cover_down).any()

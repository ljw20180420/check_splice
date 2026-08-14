import os
import pathlib

from .check import ReadStart, all_intervals
from .sam import filter_reads, parse_block


def process_locus(
    samfile: os.PathLike,
    chrom: str,
    start: int,
    end: int,
    cpcdh_file: os.PathLike,
    cover_threshold: int = 3,
):
    samfile = pathlib.Path(os.fspath(samfile))
    exp, protein, clone, rep = samfile.stem.split("_")
    introns, intron_starts, intron_ends = all_intervals(cpcdh_file, cover_threshold)
    for read in filter_reads(samfile, chrom, start, end):
        info = {
            "exp": exp,
            "protein": protein,
            "clone": clone,
            "rep": rep,
            "query_name": read.query_name,
            "is_read1": read.is_read1,
            "is_qcfail": read.is_qcfail,
            "is_duplicate": read.is_duplicate,
            "mapping_quality": read.mapping_quality,
            "is_shadow": False,
        }

        blocks = list(parse_block(read))
        info = introns(info, blocks, "connect")
        info = intron_starts(info, blocks, "cover")
        info = intron_ends(info, blocks, "cover")
        info = ReadStart.get(info, blocks, read.is_read1)

        yield info

        if exp == "rna" and protein in ["MPP8", "PPHLN1", "TASOR"]:
            info_shadow = info.copy()
            info_shadow["is_shadow"] = True
            read.is_read1 = not read.is_read1

            blocks = list(parse_block(read))
            info_shadow = introns(info_shadow, blocks, "connect")
            info_shadow = intron_starts(info_shadow, blocks, "cover")
            info_shadow = intron_ends(info_shadow, blocks, "cover")
            info_shadow = ReadStart.get(info_shadow, blocks, read.is_read1)

            yield info_shadow


def process_all(
    bam_dir: os.PathLike,
    chrom: str,
    start: int,
    end: int,
    cpcdh_file: os.PathLike,
    cover_threshold: int = 3,
):
    bam_dir = pathlib.Path(os.fspath(bam_dir))
    for samfile in os.listdir(bam_dir):
        if not samfile.endswith(".bam"):
            continue

        yield from process_locus(
            samfile=bam_dir / samfile,
            chrom=chrom,
            start=start,
            end=end,
            cpcdh_file=cpcdh_file,
            cover_threshold=cover_threshold,
        )

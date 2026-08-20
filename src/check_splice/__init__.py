import os
import pathlib

from .check import Read, all_intervals
from .sam import filter_reads, parse_block


def process_locus(
    samfile: os.PathLike,
    chrom: str,
    start: int,
    end: int,
    cpcdh_file: os.PathLike,
    cover_threshold: int,
    exon_end_extend: int,
):
    samfile = pathlib.Path(os.fspath(samfile))
    exp, protein, clone, rep = samfile.stem.split("_")
    (
        introns,
        intron_starts,
        intron_ends,
        exon_ends,
    ) = all_intervals(cpcdh_file, cover_threshold, exon_end_extend)

    for read in filter_reads(samfile, chrom, start, end):
        info = {
            "exp": exp,
            "protein": protein,
            "clone": clone,
            "rep": rep,
            "query_name": read.query_name,
            "is_forward": read.is_forward,
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
        info = exon_ends(info, blocks, "inrange_end")
        info = Read.start(info, blocks, read.is_read1)

        yield info

        if exp == "rna" and protein in ["MPP8", "PPHLN1", "TASOR"]:
            info_shadow = info.copy()
            info_shadow["is_shadow"] = True
            read.is_read1 = not read.is_read1

            blocks = list(parse_block(read))
            info_shadow = introns(info_shadow, blocks, "connect")
            info_shadow = intron_starts(info_shadow, blocks, "cover")
            info_shadow = intron_ends(info_shadow, blocks, "cover")
            info_shadow = Read.start(info_shadow, blocks, read.is_read1)

            yield info_shadow


def process_all(cfg: dict):
    bam_dir = cfg["data_dir"] / "bam"
    chrom = cfg["chrom"]
    start = cfg["start"]
    end = cfg["end"]
    cpcdh_file = cfg["data_dir"] / "result" / "cpcdh.csv"
    cover_threshold = cfg["cover_threshold"]
    exon_end_extend = cfg["exon_end_extend"]

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
            exon_end_extend=exon_end_extend,
        )

import os
import pathlib

from .check import Read, all_intervals
from .sam import filter_reads, parse_block_with_flip


def process_locus(
    bamfile: os.PathLike,
    chrom: str,
    start: int,
    end: int,
    cpcdh_file: os.PathLike,
    cover_threshold: int,
    exon_end_extend: int,
    flip: str,
):
    bamfile = pathlib.Path(os.fspath(bamfile))
    exp, protein, clone, rep = bamfile.stem.split("_")
    (
        introns,
        intron_starts,
        intron_ends,
        exon_ends,
    ) = all_intervals(cpcdh_file, cover_threshold, exon_end_extend)

    info = {
        "exp": exp,
        "protein": protein,
        "clone": clone,
        "rep": rep,
    }
    for read in filter_reads(bamfile, chrom, start, end):
        info |= {
            "query_name": read.query_name,
            "is_forward": read.is_forward,
            "is_read1": read.is_read1,
            "is_qcfail": read.is_qcfail,
            "is_duplicate": read.is_duplicate,
            "mapping_quality": read.mapping_quality,
            "is_shadow": False,
        }

        blocks = parse_block_with_flip(read, flip)
        info["blocks"] = ";".join([
            f"{block_chrom}:{block_start}:{block_end}:{block_strand}"
            for block_chrom, block_start, block_end, block_strand in blocks
        ])
        info = introns(info, blocks, "connect")
        info = intron_starts(info, blocks, "cover")
        info = intron_ends(info, blocks, "cover")
        info = exon_ends(info, blocks, "inrange_end")
        info = Read.start(info, blocks, read.is_read1)

        yield info.copy()

        if exp == "rna" and protein in ["MPP8", "PPHLN1", "TASOR"]:
            info["is_shadow"] = True
            read.is_read1 = not read.is_read1

            blocks = parse_block_with_flip(read, flip)
            info["blocks"] = ";".join([
                f"{block_chrom}:{block_start}:{block_end}:{block_strand}"
                for block_chrom, block_start, block_end, block_strand in blocks
            ])
            info = introns(info, blocks, "connect")
            info = intron_starts(info, blocks, "cover")
            info = intron_ends(info, blocks, "cover")
            info = Read.start(info, blocks, read.is_read1)

            yield info.copy()


def process_all(cfg: dict, assemble: str):
    bam_dir = cfg["data_dir"] / "bam"
    chrom = cfg[assemble]["chrom"]
    start = cfg[assemble]["start"]
    end = cfg[assemble]["end"]
    cpcdh_file = cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
    cover_threshold = cfg["cover_threshold"]
    exon_end_extend = cfg["exon_end_extend"]

    for bamfile in os.listdir(bam_dir):
        if not bamfile.endswith(".bam"):
            continue

        exp, protein, clone, rep = bamfile.removesuffix(".bam").split("_")
        if assemble == "hg19":
            if clone.startswith("mm"):
                continue
        elif assemble == "mm10":
            if not clone.startswith("mm"):
                continue
        else:
            raise ValueError("unknown assemble")

        yield from process_locus(
            bamfile=bam_dir / bamfile,
            chrom=chrom,
            start=start,
            end=end,
            cpcdh_file=cpcdh_file,
            cover_threshold=cover_threshold,
            exon_end_extend=exon_end_extend,
            flip=cfg["flip"],
        )

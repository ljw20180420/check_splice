import os
import pathlib

import pysam
import sh

from .common import ParseSamRead, filter_bam_reads
from .utils import clone2assemble, get_precursor_pos, get_sample_bam


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
                    ) = parse_sam_read.flip_read(
                        ref_block_chrom,
                        ref_block_start,
                        ref_block_end,
                        ref_block_strand,
                        query_block_start,
                        query_block_end,
                        align_string,
                        read.query_length,
                    )

                yield f"{exp},{protein},{clone},{rep},{read.query_name},{read.get_forward_sequence()}{read.is_forward},{read.is_read1},{read.is_qcfail},{read.is_duplicate},{read.mapping_quality},{ref_block_chrom},{ref_block_start},{ref_block_end},{ref_block_strand},{query_block_start},{query_block_end},{align_string}"


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
    parse_sam_read = ParseSamRead()
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
                    ref_block_chrom,
                    ref_block_start,
                    ref_block_end,
                    ref_block_strand,
                    query_block_start,
                    query_block_end,
                    align_string,
                ) in parse_sam_read.parse_read(read):
                    cover_up = (df_se["pos"] + cfg["cover_threshold"]).between(
                        ref_block_start, ref_block_end
                    )
                    cover_down = (df_se["pos"] - cfg["cover_threshold"]).between(
                        ref_block_start, ref_block_end
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

                if "N" in read.cigarstring or read.has_tag("SA"):
                    outfile.write(read)

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

import os
import re

import pandas as pd
import pyBigWig
import pysam
import sh


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


def get_precursor_pos(cfg: dict) -> pd.DataFrame:
    df_se = (
        pd
        .read_csv(cfg["data_dir"] / "result" / "cpcdh.csv")
        .query("type == 'exon'")
        .melt(
            id_vars=["chrom", "name"],
            value_vars=["start", "end"],
            var_name="se",
            value_name="pos",
        )
    )

    return df_se


def modify_path(bam_file: os.PathLike, RS: str) -> os.PathLike:
    return bam_file.with_name("precursor") / bam_file.with_suffix(f".{RS}.bam").name


def filter_precursor_reads(cfg: dict, bam_file: os.PathLike, RS: str) -> None:
    df_se = get_precursor_pos(cfg)
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        with pysam.AlignmentFile(
            modify_path(bam_file, RS), "wb", template=infile
        ) as outfile:
            for read in infile.fetch(cfg["chrom"], cfg["start"], cfg["end"]):
                if read.is_secondary:
                    continue
                if not read.is_mapped:
                    continue
                if read.is_supplementary:
                    continue
                if ("1" in RS) != read.is_read1:
                    continue
                if ("f" in RS) != read.is_forward:
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


def filter_precursor_reads_all(cfg: dict) -> None:
    (cfg["data_dir"] / "bam" / "precursor").mkdir(parents=True, exist_ok=True)
    samtools = sh.Command("samtools")
    for bam_file in os.listdir(cfg["data_dir"] / "bam"):
        if not bam_file.endswith(".bam"):
            continue

        bam_file = cfg["data_dir"] / "bam" / bam_file
        for RS in ["1r", "2f", "1f", "2r"]:
            filter_precursor_reads(cfg, bam_file, RS)

        for S, RS1, RS2 in zip(["f", "r"], ["1r", "1f"], ["2f", "2r"]):
            samtools(
                "merge",
                "-f",
                "-o",
                os.fspath(modify_path(bam_file, S)),
                os.fspath(modify_path(bam_file, RS1)),
                os.fspath(modify_path(bam_file, RS2)),
            )


def merge_bam_by_exp_protein_treat(cfg: dict) -> None:
    df_total_count = (
        pd
        .read_csv(cfg["data_dir"] / "result" / "total_count.csv", header=0)
        .assign(
            wt=lambda df: df["clone"].map(
                lambda ele: "control" if ele.startswith("WT") else "treat"
            )
        )
        .groupby(["exp", "protein", "wt"])
        .agg(total_count=pd.NamedAgg("total_count", "sum"))
    ).reset_index()
    (cfg["data_dir"] / "bam" / "precursor" / "merge").mkdir(exist_ok=True, parents=True)
    samtools = sh.Command("samtools")
    bamCoverage = sh.Command("bamCoverage")
    for exp in ["total", "rna", "pro", "clip"]:
        for protein in ["WT", "NP220", "MPP8", "PPHLN1", "TASOR"]:
            for wt in [True, False]:
                if wt:
                    treat = "control"
                    total_counts = df_total_count.query(
                        "exp == @exp and protein == @protein and wt == 'control'"
                    )["total_count"]
                else:
                    if exp != "clip":
                        treat = "delta"
                    else:
                        treat = "tag"
                    total_counts = df_total_count.query(
                        "exp == @exp and protein == @protein and wt == 'treat'"
                    )["total_count"]
                assert len(total_counts) < 2, "more than one total count found"
                total_count = total_counts.sum().item()

                for strand in ["f", "r"]:
                    bam_files = []
                    for bam_file in os.listdir(cfg["data_dir"] / "bam" / "precursor"):
                        if not bam_file.endswith(f".{strand}.bam"):
                            continue
                        exp_, protein_, clone_, _ = bam_file.split("_", 3)
                        if exp_ != exp or protein_ != protein:
                            continue
                        if wt != clone_.startswith("WT"):
                            continue

                        bam_file = cfg["data_dir"] / "bam" / "precursor" / bam_file
                        bam_files.append(os.fspath(bam_file))

                    if not bam_files:
                        continue

                    merge_bam = (
                        cfg["data_dir"]
                        / "bam"
                        / "precursor"
                        / "merge"
                        / f"{exp}_{protein}_{treat}.{strand}.bam"
                    )
                    samtools(
                        "merge",
                        "-f",
                        "-o",
                        os.fspath(merge_bam),
                        *bam_files,
                    )
                    samtools("index", os.fspath(merge_bam))

                    with pysam.AlignmentFile(merge_bam, "rb") as bam:
                        mapped_count = bam.mapped
                    if mapped_count > 0:
                        bamCoverage(
                            "--bam",
                            os.fspath(merge_bam),
                            "-o",
                            os.fspath(merge_bam.with_suffix(".bw")),
                            "-r",
                            f"{cfg['chrom']}:{cfg['start']}:{cfg['end']}",
                            "--binSize",
                            10,
                            "--scaleFactor",
                            f"{1_000_000 / total_count}",
                        )
                    else:
                        mid = (cfg["start"] + cfg["end"]) // 2
                        with pyBigWig.open(
                            os.fspath(merge_bam.with_suffix(".bw")), "w"
                        ) as bw:
                            bw.addHeader([("chr5", 180915260)])
                            bw.addEntries("chr5", [mid], values=[0.0], span=1)

                    with pysam.AlignmentFile(merge_bam, "rb") as ib:
                        with pysam.AlignmentFile(
                            merge_bam.with_suffix(".unify.bam"), "wb", template=ib
                        ) as ob:
                            for read in ib.fetch():
                                read.is_forward = strand == "f"
                                ob.write(read)

                    samtools("index", os.fspath(merge_bam.with_suffix(".unify.bam")))

import os
import pathlib
import shutil
import subprocess

import numpy as np
import pandas as pd
import pyBigWig
import pysam
import sh


def get_sample_bam(cfg: dict) -> pd.DataFrame:
    exps = []
    proteins = []
    clones = []
    reps = []
    files = []
    for bamfile in os.listdir(cfg["data_dir"] / "bam"):
        if not bamfile.endswith(".bam"):
            continue

        exp, protein, clone, rep = bamfile.removesuffix(".bam").split("_")
        exps.append(exp)
        proteins.append(protein)
        clones.append(clone)
        reps.append(rep)
        files.append(os.fspath(cfg["data_dir"] / "bam" / bamfile))

    return pd.DataFrame({
        "exp": exps,
        "protein": proteins,
        "clone": clones,
        "rep": reps,
        "file": files,
    })


def clone2assemble(clone: str) -> str:
    if clone.startswith("mm"):
        return "mm10"
    return "hg19"


def clone2treat(clone: str) -> str:
    if clone.startswith("mmWT"):
        return "mmcontrol"
    if clone.startswith("WT"):
        return "control"
    treat = "treat"
    if clone.startswith("mm"):
        treat = f"mm{treat}"

    return treat


def blocks_string2tuple(blocks: str):
    for block in blocks.split(";"):
        chrom, start, end, strand = block.split(":")
        yield chrom, int(start), int(end), strand


class SelectTotalCount:
    def __init__(self, cfg: dict) -> None:
        self.df_total = (
            pd
            .read_csv(cfg["data_dir"] / "result" / "total_count.csv", header=0)
            .assign(treat=lambda df: df["clone"].map(clone2treat))
            .groupby(by=["exp", "protein", "treat"], as_index=False)["total_count"]
            .sum()
        )

    def __call__(self, exp: str, protein: str, treat: str) -> int:
        return self.df_total.query(
            "exp == @exp and protein == @protein and treat == @treat"
        )["total_count"].item()


def get_bw(
    bamfile: os.PathLike,
    total_count: int,
    bin_size: int,
    chrom: str,
    start: int,
    end: int,
    chrom_size,
) -> None:
    bamfile = pathlib.Path(bamfile)
    bamCoverage = sh.Command("bamCoverage")
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        mapped_count = bam.mapped
    if mapped_count > 0:
        bamCoverage(
            "--bam",
            os.fspath(bamfile),
            "-o",
            os.fspath(bamfile.with_suffix(".bw")),
            "-r",
            f"{chrom}:{start}:{end}",
            "--binSize",
            bin_size,
            "--scaleFactor",
            1_000_000 / total_count,
        )
    else:
        mid = (start + end) // 2
        with pyBigWig.open(os.fspath(bamfile.with_suffix(".bw")), "w") as bw:
            bw.addHeader([(chrom, chrom_size)])
            bw.addEntries(chrom, [mid], values=[0.0], span=1)


def bw_merge_adjacent_intervals_with_identical_values(
    starts: np.ndarray, ends: np.ndarray, values: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    change_indices = np.where(np.diff(values) != 0)[0] + 1

    ends = np.concatenate((starts[change_indices], [ends[-1]]))
    starts = starts[np.concatenate(([0], change_indices))]
    values = values[np.concatenate(([0], change_indices))]

    return starts, ends, values


def get_precursor_pos(cfg: dict, assemble: str) -> pd.DataFrame:
    df_se = (
        pd
        .read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv")
        .query("type == 'exon'")
        .melt(
            id_vars=["chrom", "name"],
            value_vars=["start", "end"],
            var_name="se",
            value_name="pos",
        )
    )

    return df_se

import os
import pathlib
import shutil
import subprocess

import numpy as np
import pandas as pd
import pyarrow as pa
import pyBigWig
import pysam
import sh
from pyarrow import ipc


def get_total_count(cfg: dict) -> None:
    exps = []
    proteins = []
    clones = []
    reps = []
    total_counts = []
    for bamfile in os.listdir(cfg["data_dir"] / "bam"):
        if not bamfile.endswith(".bam"):
            continue

        exp, protein, clone, rep = bamfile.removesuffix(".bam").split("_")
        bamfile = cfg["data_dir"] / "bam" / bamfile
        with pysam.AlignmentFile(os.fspath(bamfile)) as bam:
            total_count = sum(
                1
                for read in bam
                if not read.is_secondary
                and read.is_mapped
                and not read.is_supplementary
            )

        print(bamfile, total_count)

        exps.append(exp)
        proteins.append(protein)
        clones.append(clone)
        reps.append(rep)
        total_counts.append(total_count)

    pd.DataFrame({
        "exp": exps,
        "protein": proteins,
        "clone": clones,
        "rep": reps,
        "total_count": total_counts,
    }).to_csv(cfg["data_dir"] / "result" / "total_count.csv", index=False)


def select_total_count(cfg: dict, exp: str, protein: str, treat: str) -> int:
    df_total = pd.read_csv(cfg["data_dir"] / "result" / "total_count.csv", header=0)
    df_total = (
        df_total
        .assign(treat=lambda df: get_treat(df))
        .groupby(["exp", "protein", "treat"])["total_count"]
        .sum()
        .reset_index()
    )

    total_count = df_total.query(
        "exp == @exp and protein == @protein and treat == @treat"
    )["total_count"].item()

    return total_count


def jsonl2feather(jsonl_file: os.PathLike, feather_file: os.PathLike):
    writer = None
    # Stream JSON in chunks of 10,000 rows
    for chunk in pd.read_json(jsonl_file, chunksize=100000, lines=True):
        batch = pa.RecordBatch.from_pandas(chunk)
        if writer is None:
            writer = ipc.RecordBatchFileWriter(feather_file, batch.schema)
        writer.write_batch(batch)

    if writer is not None:
        writer.close()

    df = pd.read_feather(feather_file)
    df.to_feather(feather_file)


def pair_to_hic(
    pair_file: os.PathLike, resolutions: list[int], chrom_sizes: os.PathLike
) -> None:
    hic_file = pair_file.with_suffix(".hic")
    subprocess.run(
        args=[
            "hictk",
            "load",
            "--format",
            "4dn",
            "--bin-size",
            f"{resolutions[0]}",
            "--chrom-sizes",
            os.fspath(chrom_sizes),
            "--force",
            os.fspath(pair_file),
            os.fspath(hic_file),
        ],
        check=False,
    )

    subprocess.run(
        args=["hictk", "zoomify", "--resolutions"]
        + [f"{resolution}" for resolution in resolutions]
        + [
            "--force",
            os.fspath(hic_file),
            f"{os.fspath(hic_file.with_suffix('.m.hic'))}",
        ],
        check=False,
    )

    shutil.move(f"{os.fspath(hic_file.with_suffix('.m.hic'))}", os.fspath(hic_file))

    subprocess.run(args=["hictk", "balance", "scale", os.fspath(hic_file)], check=False)


def get_treat(df: pd.DataFrame) -> pd.Series:
    treat = df["clone"].map(
        lambda ele: (
            "mmcontrol"
            if ele.startswith("mmWT")
            else "mmdelta"
            if ele.startswith("mm")
            else "control"
            if ele.startswith("WT")
            else "delta"
        )
    )
    treat = treat.where((df["exp"] != "clip") | (treat.str.contains("control")), "tag")
    treat = treat.where((treat != "tag") | (~df["clone"].str.startswith("mm")), "mmtag")

    return treat


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

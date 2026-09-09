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


def prepare_gene_bed12(cfg: dict) -> None:
    gtfToGenePred = sh.Command("gtfToGenePred")
    gtfToGenePred(
        "-genePredExt",
        os.fspath(cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gtf"),
        os.fspath(cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gp"),
    )
    df_gp = pd.read_csv(
        cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gp", sep="\t", header=None
    )
    df_gp[[11] + list(range(1, 15))].to_csv(
        cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gp",
        sep="\t",
        header=False,
        index=False,
    )
    genePredToBigGenePred = sh.Command("genePredToBigGenePred")
    genePredToBigGenePred(
        os.fspath(cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gp"),
        os.fspath(cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.bgp"),
    )
    df_bgp = (
        pd
        .read_csv(
            cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.bgp", sep="\t", header=None
        )[list(range(12))]
        .rename(
            columns={
                0: "chrom",
                1: "chromStart",
                2: "chromEnd",
                3: "name",
                4: "score",
                5: "strand",
                6: "thickStart",
                7: "thickEnd",
                8: "itemRgb",
                9: "blockCount",
                10: "blockSizes",
                11: "blockStarts",
            }
        )
        .sort_values(by=["chrom", "chromStart"], ignore_index=True)
        .query("not name.str.startswith('PCDHA') or blockCount == 4")
        .query(
            "not name.str.startswith('PCDHB') or blockCount == 1 or name == 'PCDHB9'"
        )
        .query("not name.str.startswith('PCDHB@')")
        .query("not name.str.startswith('PCDHG') or blockCount == 4")
        .query("name != 'PCDHA1' or blockSizes.str.startswith('2545')")
        .query("name != 'PCDHA6' or blockSizes.str.startswith('2526')")
        .query("name != 'PCDHA10' or blockSizes.str.startswith('2540')")
        .query("name != 'PCDHGA11' or blockSizes.str.startswith('2610')")
        .query("name != 'PCDHGC3' or blockSizes.str.startswith('2581')")
        .query(
            "name != 'LOC112267934' and name != 'LOC101926905' and name != 'LOC100419552' and name != 'SLC25A2' and name != 'TAF7' and name != 'RN7SL68P'"
        )
        .reset_index(drop=True)
    )

    df_bgp.to_csv(
        cfg["data_dir"] / "result" / "hg19.12.bed", sep="\t", header=False, index=False
    )
    (cfg["data_dir"] / "result" / "hg19.12.bed.bgz").unlink(missing_ok=True)
    (cfg["data_dir"] / "result" / "hg19.12.bed.bgz.tbi").unlink(missing_ok=True)


def get_treat(df: pd.DataFrame) -> pd.Series:
    treat = df["clone"].map(lambda ele: "control" if ele.startswith("WT") else "delta")
    treat = treat.where((df["exp"] != "clip") | (treat == "control"), "tag")

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


def get_total_count(cfg: dict, exp: str, protein: str, treat: str) -> int:
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

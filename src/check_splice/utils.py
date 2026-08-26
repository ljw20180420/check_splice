import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import pyarrow as pa
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


def write_hic(
    df: pd.DataFrame,
    hic_file: os.PathLike,
    resolutions: list[int],
    chrom_sizes: os.PathLike,
) -> None:
    pair_file = hic_file.with_suffix(".pairs")
    with open(pair_file, "w") as fd:
        fd.write("## pairs format v1.0\n")
        df.rename(
            columns={
                "pos1": "pos1_old",
                "pos2": "pos2_old",
            }
        ).assign(
            pos1=lambda df: np.minimum(df["pos1_old"], df["pos2_old"]),
            pos2=lambda df: np.maximum(df["pos1_old"], df["pos2_old"]),
        )[["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]].to_csv(
            fd, sep="\t", header=False, index=False
        )

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
    df = (
        pd
        .read_csv(
            "/home/ljw/sdb1/ucsc/hubs/myHub/lhg19/lhg19.bgp", sep="\t", header=None
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
        .query("not name.str.startswith('PCDHA') or blockCount == 4")
        .query(
            "not name.str.startswith('PCDHB') or blockCount == 1 or name == 'PCDHB9'"
        )
        .query("not name.str.startswith('PCDHG') or blockCount == 4")
        .query("name != 'PCDHA1' or blockSizes.str.startswith('2545')")
        .query("name != 'PCDHA6' or blockSizes.str.startswith('2526')")
        .query("name != 'PCDHA10' or blockSizes.str.startswith('2540')")
        .query("name != 'PCDHGA11' or blockSizes.str.startswith('2610')")
        .query("name != 'PCDHGC3' or blockSizes.str.startswith('2581')")
        .reset_index(drop=True)
    )

    df.to_csv(
        cfg["data_dir"] / "result" / "hg19.12.bed", sep="\t", header=False, index=False
    )

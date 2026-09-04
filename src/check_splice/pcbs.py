import os

import pandas as pd


def get_pCBS(cfg: dict) -> None:
    shift_file = "pCBS_shift.csv"
    cpcdh_file = cfg["data_dir"] / "result" / "cpcdh.csv"

    df_shift = pd.read_csv(shift_file, header=0).assign(
        name=lambda df: df["gene"].str.upper()
    )
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    df = (
        pd
        .merge(left=df_shift, right=df_cpcdh, how="inner", on="name", validate="1:1")
        .astype({"CDS_start": float, "CDS_end": float})
        .astype({"CDS_start": int, "CDS_end": int})
    )
    df = df.assign(
        end=lambda df: df["CDS_start"] + df["shift"] + 2,
        start=lambda df: df["end"] - 27,
    )[["chrom", "start", "end", "name"]].assign(score=".", strand="+")

    df.to_csv(
        cfg["data_dir"] / "result" / "pCBS.bed", sep="\t", header=False, index=False
    )
    (cfg["data_dir"] / "result" / "pCBS.bed.bgz").unlink(missing_ok=True)
    (cfg["data_dir"] / "result" / "pCBS.bed.bgz.tbi").unlink(missing_ok=True)

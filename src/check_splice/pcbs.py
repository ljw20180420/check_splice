import os

import pandas as pd


def get_pCBS(shift_file: os.PathLike, cpcdh_file: os.PathLike) -> pd.DataFrame:
    df_shift = pd.read_csv(shift_file, header=0).assign(
        name=lambda df: df["gene"].str.upper()
    )
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    df = pd.merge(left=df_shift, right=df_cpcdh, how="inner", on="name", validate="1:1")
    df = df.assign(
        start=lambda df: df["CDS_start"] + df["shift"],
        end=lambda df: df["start"] + 42,
    )[["chrom", "start", "end", "name"]].assign(score=".", strand="+")

    return df

import os

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .draw import _heatmap


def splice(result_file: os.PathLike, cpcdh_file: os.PathLike) -> pd.DataFrame:
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    intron_names = df_cpcdh.query("type='intron'")["name"].to_list()

    df = pd.read_json(result_file, lines=True)
    columns = [
        column
        for intro_name in intron_names
        for column in (
            f"connect_{intro_name}",
            f"cover_{intro_name}_start",
            f"cover_{intro_name}_end",
        )
    ]
    df = (
        df
        .groupby(["exp", "protein", "clone", "rep", "query_name"])
        .agg(**{
            column: pd.NamedAgg(column=column, aggfunc="any") for column in columns
        })
        .reset_index()
    )

    df = (
        df
        .groupby([
            "exp",
            "clone",
        ])
        .agg(**{
            column: pd.NamedAgg(column=column, aggfunc="sum") for column in columns
        })
        .reset_index()
    )

    df = (
        df
        .melt(
            id_vars=["exp", "clone"],
            value_vars=columns,
            var_name=["opt_intron"],
            value_name="count",
        )
        .assign(
            intron=lambda df: df["opt_intron"].str.rsplit("_", n=1, expand=True)[1],
            opt=lambda df: df["opt_intron"].str.rsplit("_", n=1, expand=True)[0],
        )
        .pivot_table(values="count", index=["exp", "clone", "intron"], columns="opt")
        .reset_index()
    )

    return df


def around_tss(
    result_file: os.PathLike, cpcdh_file: os.PathLike, tss_extend: int
) -> tuple[Figure, Axes]:
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    tsses = df_cpcdh.query("type='exon' and name.str.startswith('PCDH')")[
        ["start", "name"]
    ]

    df = pd.read_json(result_file, lines=True)
    read_starts = (
        df
        .query("not is_shadow and is_read1 and strand == '+'")
        .reset_index(drop=True)["read_start"]
        .value_counts()
    )
    arr = np.stack(
        [
            read_starts
            .reindex(range(tss - tss_extend, tss + tss_extend + 1), fill_value=0)
            .sort_index()
            .to_numpy()
            for tss in tsses["start"]
        ],
        axis=0,
    )

    fig, ax = _heatmap(mat=arr, extent=[-tss_extend, tss_extend, 0, len(tsses) - 1])
    ax.set_xlabel("position")
    ax.set_ylabel("exon")
    ax.set_yticklabels(labels=tsses["name"].tolist())

    return fig, ax

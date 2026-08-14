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
        .reset_index(drop=True)[["exp", "clone", "read_start"]]
        .value_counts()
        .reset_index()
    )

    exps = read_starts["exp"].drop_duplicates().tolist()
    clones = read_starts["clone"].drop_duplicates().tolist()

    df_arounds = []
    for tss, name in zip(tsses["start"], tsses["name"]):
        df_arounds.append(
            read_starts
            .assign(relative=lambda df, tss=tss: df["read_start"] - tss)
            .query("relative >= -@tss_extend and relative <= @tss_extend")
            .assign(exon=name)
        )
    df_around = (
        pd
        .concat(df_arounds, ignore_index=True)
        .pivot_table(values="count", index=["exp", "clone", "exon"], columns="relative")
        .reindex(
            index=pd.MultiIndex.from_product(
                [
                    exps,
                    clones,
                    tsses["name"].tolist(),
                ],
                names=["exp", "clone", "exon"],
            ),
            columns=range(-tss_extend, tss_extend + 1),
            fill_value=0,
        )
        .reset_index()
    )

    fig, ax = _heatmap(
        mat=df_around.to_numpy(), extent=[-tss_extend, tss_extend, 0, len(tsses) - 1]
    )
    ax.set_xlabel("position")
    ax.set_ylabel("exon")
    ax.set_yticklabels(labels=df.index.map(lambda tup: "_".join(str(e) for e in tup)))

    return fig, ax

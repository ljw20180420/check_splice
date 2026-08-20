import matplotlib
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_text,
    geom_raster,
    ggplot,
    scale_fill_gradient,
    scale_y_discrete,
    theme,
)

matplotlib.use("agg")


def swap_elements(vec: list, a: str, b: str) -> list:
    i, j = vec.index(a), vec.index(b)
    vec[j], vec[i] = vec[i], vec[j]
    return vec


def splice(cfg: dict) -> pd.DataFrame:
    cpcdh_file = cfg["data_dir"] / "result" / "cpcdh.csv"
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    intron_names = df_cpcdh.query("type=='intron'")["name"].to_list()

    read_file = cfg["data_dir"] / "result" / "reads.feather"
    df = pd.read_feather(read_file)

    mask = np.any(
        np.array([
            df.columns.str.endswith(f".{intro_name}") for intro_name in intron_names
        ]),
        axis=0,
    )
    columns = df.columns[mask].to_list()

    df = (
        df
        .groupby(["exp", "protein", "clone", "rep", "query_name"])
        .agg(**{
            column: pd.NamedAgg(column=column, aggfunc="any") for column in columns
        })
        .copy()
        .reset_index()
    )

    df_total = pd.read_csv(cfg["data_dir"] / "result" / "total_count.csv", header=0)
    df_total = (
        df_total
        .assign(is_WT=lambda df: df["clone"].str.startswith("WT"))
        .groupby(["exp", "protein", "is_WT"])["total_count"]
        .sum()
        .reset_index()
    )

    df = (
        df
        .assign(is_WT=lambda df: df["clone"].str.startswith("WT"))
        .groupby(["exp", "protein", "is_WT"])
        .agg(**({
            column: pd.NamedAgg(column=column, aggfunc="sum") for column in columns
        }))
        .copy()
        .reset_index()
        .merge(
            df_total, how="left", on=["exp", "protein", "is_WT"], validate="many_to_one"
        )
    )

    df = (
        df
        .melt(
            id_vars=["exp", "protein", "is_WT", "total_count"],
            value_vars=columns,
            var_name="opt_intron",
            value_name="count",
        )
        .assign(
            intron=lambda df: df["opt_intron"].str.rsplit(".", n=1, expand=True)[1],
            opt=lambda df: df["opt_intron"].str.rsplit(".", n=1, expand=True)[0],
        )
        .pivot_table(
            values="count",
            index=["exp", "protein", "is_WT", "total_count", "intron"],
            columns="opt",
        )
        .reset_index()
    )

    df = (
        pd
        .merge(
            df,
            df_cpcdh[["chrom", "start", "end", "name"]],
            how="left",
            left_on="intron",
            right_on="name",
            validate="many_to_one",
        )
        .sort_values(by=["exp", "protein", "is_WT", "start"], ignore_index=True)
        .drop(columns="intron")
    )

    df = df[swap_elements(df.columns.to_list(), "cover.start", "cover.end")]

    return df


def around_tss(
    cfg: dict,
) -> None:
    cpcdh_file = cfg["data_dir"] / "result" / "cpcdh.csv"
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    tsses = df_cpcdh.query("type=='exon' and name.str.startswith('PCDH')")[
        ["start", "name"]
    ]

    result_file = cfg["data_dir"] / "result" / "reads.feather"
    df = pd.read_feather(result_file).assign(
        is_WT=lambda df: df["clone"].str.startswith("WT"),
        exp_protein_wt=lambda df: (
            df["exp"]
            + "_"
            + df["protein"]
            + "_"
            + df["is_WT"].map({True: "wt", False: "non"})
        ),
    )
    exp_protein_wts = df["exp_protein_wt"].drop_duplicates().to_list()
    read_starts = (
        df
        .query("not is_shadow and is_read1 and is_forward")
        .reset_index(drop=True)[["exp_protein_wt", "read_start"]]
        .value_counts()
        .reset_index()
    )

    tss_extend = cfg["tss_extend"]
    df_arounds = []
    for tss, name in zip(tsses["start"], tsses["name"]):
        df_arounds.append(
            read_starts
            .assign(relative=lambda df, tss=tss: df["read_start"] - tss)
            .query("relative >= -@tss_extend and relative <= @tss_extend")
            .assign(exon=name)
        )
    df_around = pd.concat(df_arounds, ignore_index=True)
    df_around["exon"] = pd.Categorical(
        df_around["exon"], categories=tsses["name"].to_list(), ordered=True
    )

    df_around_agg = (
        df_around
        .groupby(["exon", "relative"])
        .agg(count=pd.NamedAgg(column="count", aggfunc="sum"))
        .reindex(
            index=pd.MultiIndex.from_product(
                [
                    tsses["name"].to_list(),
                    list(range(-tss_extend, tss_extend + 1)),
                ],
                names=["exon", "relative"],
            ),
            fill_value=0,
        )
        .reset_index()
    )
    (
        ggplot(df_around_agg, mapping=aes(x="relative", y="exon", fill="count"))
        + geom_raster()
        + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
        + scale_y_discrete(limits=tsses["name"].to_list()[::-1])
        + theme(axis_text_x=element_text(angle=90, ma="right"))
    ).save(cfg["data_dir"] / "result" / "agg_around_tss.png")

    for exp_protein_wt in exp_protein_wts:
        df_slice = (
            df_around
            .query("exp_protein_wt == @exp_protein_wt")
            .reset_index(drop=True)[["exon", "relative", "count"]]
            .set_index(["exon", "relative"])
            .reindex(
                index=pd.MultiIndex.from_product(
                    [
                        tsses["name"].to_list(),
                        list(range(-tss_extend, tss_extend + 1)),
                    ],
                    names=["exon", "relative"],
                ),
                fill_value=0,
            )
            .reset_index()
        )
        (
            ggplot(df_slice, mapping=aes(x="relative", y="exon", fill="count"))
            + geom_raster()
            + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
            + scale_y_discrete(limits=tsses["name"].to_list()[::-1])
            + theme(axis_text_x=element_text(angle=90, ma="right"))
        ).save(cfg["data_dir"] / "result" / f"{exp_protein_wt}_around_tss.png")

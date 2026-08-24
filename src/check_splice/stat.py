from collections.abc import Iterable

import matplotlib
import numpy as np
import pandas as pd
import pypdf

from .draw import around_heatmap
from .utils import write_hic

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


def around(
    cfg: dict,
    centers: Iterable[int],
    center_names: Iterable[str],
    center_axis_name: str,
    extend: int,
    targets: list[str],
):
    result_file = cfg["data_dir"] / "result" / "reads.feather"
    df = pd.read_feather(result_file).assign(
        exp_protein_wt=lambda df: (
            df["exp"]
            + "_"
            + df["protein"]
            + "_"
            + df["clone"].map(lambda ele: "wt" if ele.startswith("WT") else "non")
        ),
    )

    df_arounds = []
    for center, center_name, target in zip(centers, center_names, targets):
        read_starts = (
            df
            .query("not is_shadow and is_read1 and is_forward")
            .reset_index(drop=True)[["exp_protein_wt", target]]
            .value_counts()
            .reset_index()
        )

        df_arounds.append(
            read_starts
            .assign(
                relative=lambda df, center=center, target=target: df[target] - center
            )
            .query("relative >= -@extend and relative <= @extend")
            .assign(**{center_axis_name: center_name})
        )
    df_around = pd.concat(df_arounds, ignore_index=True).assign(**{
        center_axis_name: lambda df: pd.Categorical(
            df[center_axis_name], categories=center_names, ordered=True
        )
    })

    df_around_agg = (
        df_around
        .groupby([center_axis_name, "relative"])
        .agg(count=pd.NamedAgg(column="count", aggfunc="sum"))
        .reindex(
            index=pd.MultiIndex.from_product(
                [
                    center_names,
                    list(range(-extend, extend + 1)),
                ],
                names=[center_axis_name, "relative"],
            ),
            fill_value=0,
        )
        .reset_index()
    )

    yield df_around_agg, "agg"

    df_total = (
        pd
        .read_csv(cfg["data_dir"] / "result" / "total_count.csv", header=0)
        .assign(
            exp_protein_wt=lambda df: (
                df["exp"]
                + "_"
                + df["protein"]
                + "_"
                + df["clone"].map(lambda ele: "wt" if ele.startswith("WT") else "non")
            ),
        )
        .groupby("exp_protein_wt")["total_count"]
        .sum()
        .reset_index()
    )

    for exp_protein_wt, total_count in zip(
        df_total["exp_protein_wt"], df_total["total_count"]
    ):
        df_slice = (
            df_around
            .query("exp_protein_wt == @exp_protein_wt")
            .reset_index(drop=True)[[center_axis_name, "relative", "count"]]
            .set_index([center_axis_name, "relative"])
            .reindex(
                index=pd.MultiIndex.from_product(
                    [
                        center_names,
                        list(range(-extend, extend + 1)),
                    ],
                    names=[center_axis_name, "relative"],
                ),
                fill_value=0,
            )
            .reset_index()
        )

        yield df_slice, exp_protein_wt

        df_slice = df_slice.assign(
            count=lambda df, total_count=total_count: (
                df["count"] / total_count * 1000_000
            )
        )

        yield df_slice, f"n_{exp_protein_wt}"


def read_start_around_exon_start(cfg: dict) -> None:
    cpcdh_file = cfg["data_dir"] / "result" / "cpcdh.csv"
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    tsses = df_cpcdh.query("type=='exon' and name.str.startswith('PCDH')")[
        ["start", "name"]
    ].reset_index(drop=True)

    centers = tsses["start"]
    center_names = tsses["name"]
    center_axis_name = "exon_start"
    extend = cfg["tss_extend"]
    targets = ["read_start"] * len(center_names)
    target_axis_name = "read_start"

    pdf_files = []
    with pypdf.PdfWriter() as pdf_writer:
        for df, slice in around(
            cfg, centers, center_names, center_axis_name, extend, targets
        ):
            pdf_file = around_heatmap(
                cfg, df, center_names, center_axis_name, target_axis_name, slice
            )
            pdf_writer.append(pdf_file)
            pdf_files.append(pdf_file)

        pdf_writer.write(
            cfg["data_dir"]
            / "result"
            / f"{target_axis_name}_around_{center_axis_name}.pdf"
        )

    for pdf_file in pdf_files:
        pdf_file.unlink()


def inrange_end_around_exon_end(cfg: dict) -> None:
    cpcdh_file = cfg["data_dir"] / "result" / "cpcdh.csv"
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    teses = df_cpcdh.query("type=='exon' and name.str.startswith('PCDH')")[
        ["end", "name"]
    ].reset_index(drop=True)

    centers = teses["end"]
    center_names = teses["name"]
    center_axis_name = "exon_end"
    extend = cfg["exon_end_extend"]
    targets = [f"inrange_end.end.{center_name}" for center_name in center_names]
    target_axis_name = "inrange_end"

    pdf_files = []
    with pypdf.PdfWriter() as pdf_writer:
        for df, slice in around(
            cfg, centers, center_names, center_axis_name, extend, targets
        ):
            pdf_file = around_heatmap(
                cfg, df, center_names, center_axis_name, target_axis_name, slice
            )
            pdf_writer.append(pdf_file)
            pdf_files.append(pdf_file)

        pdf_writer.write(
            cfg["data_dir"]
            / "result"
            / f"{target_axis_name}_around_{center_axis_name}.pdf"
        )

    for pdf_file in pdf_files:
        pdf_file.unlink()


def hic_4dn(cfg: dict) -> None:
    result_file = cfg["data_dir"] / "result" / "reads.feather"
    df = pd.read_feather(result_file).assign(
        exp_protein_wt=lambda df: (
            df["exp"]
            + "_"
            + df["protein"]
            + "_"
            + df["clone"].map(lambda ele: "wt" if ele.startswith("WT") else "non")
        ),
    )

    def splice_pair(ele: str) -> list:
        blocks = ele.split(";")
        pairs = []
        for i in range(len(blocks) - 1):
            block_chrom, block_start, block_end, block_strand = blocks[i].split(":")
            next_block_chrom, next_block_start, next_block_end, next_block_strand = (
                blocks[i + 1].split(":")
            )
            pairs.append(
                f"{block_chrom}:{int(block_end) + 1}:{block_strand}:{next_block_chrom}:{int(next_block_start) + 1}:{next_block_strand}"
            )
        return pairs

    df = (
        df
        .query("not is_shadow and blocks.str.contains(';')")
        .reset_index(drop=True)[["exp_protein_wt", "query_name", "blocks"]]
        .assign(blocks=lambda df: df["blocks"].map(splice_pair))
        .explode("blocks", ignore_index=True)
    )

    df = pd.concat(
        [
            df[["exp_protein_wt", "query_name"]].rename(
                columns={"query_name": "readID"}
            ),
            df["blocks"]
            .str.split(":", expand=True)
            .rename(
                columns={
                    0: "chrom1",
                    1: "pos1",
                    2: "strand1",
                    3: "chrom2",
                    4: "pos2",
                    5: "strand2",
                }
            )[["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]],
        ],
        axis=1,
    )

    exp_protein_wts = df["exp_protein_wt"].drop_duplicates().to_list()
    (cfg["data_dir"] / "result" / "hic").mkdir(exist_ok=True, parents=True)

    write_hic(
        df=df.query("strand1 == '+' and strand2 == '+'").drop(columns="exp_protein_wt"),
        hic_file=cfg["data_dir"] / "result" / "hic" / "ff.hic",
        resolutions=[1, 10, 100, 1000],
        chrom_sizes="chr5.chrom.sizes",
    )

    for exp_protein_wt in exp_protein_wts:
        write_hic(
            df=df.query(
                "exp_protein_wt == @exp_protein_wt and strand1 == '+' and strand2 == '+'"
            ).drop(columns="exp_protein_wt"),
            hic_file=cfg["data_dir"] / "result" / "hic" / f"{exp_protein_wt}_ff.hic",
            resolutions=[1, 10, 100, 1000],
            chrom_sizes="chr5.chrom.sizes",
        )

    write_hic(
        df=df.query("strand1 == '-' and strand2 == '-'").drop(columns="exp_protein_wt"),
        hic_file=cfg["data_dir"] / "result" / "hic" / "rr.hic",
        resolutions=[1, 10, 100, 1000],
        chrom_sizes="chr5.chrom.sizes",
    )

    for exp_protein_wt in exp_protein_wts:
        write_hic(
            df=df.query(
                "exp_protein_wt == @exp_protein_wt and strand1 == '-' and strand2 == '-'"
            ).drop(columns="exp_protein_wt"),
            hic_file=cfg["data_dir"] / "result" / "hic" / f"{exp_protein_wt}_rr.hic",
            resolutions=[1, 10, 100, 1000],
            chrom_sizes="chr5.chrom.sizes",
        )

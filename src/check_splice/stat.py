import os
import pathlib

import matplotlib
import pandas as pd
import pypdf
import pysam
from plotnine import (
    aes,
    element_text,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_fill_gradient,
    scale_y_discrete,
    theme,
)

from . import check
from .common import get_cpcdh_intron, is_mapped_primary_first
from .utils import (
    SelectTotalCount,
    blocks_string2tuple,
    clone2assemble,
    clone2treat,
    get_merge_bam,
    treat2assemble,
)

matplotlib.use("agg")


def splice(cfg: dict, assemble: str) -> None:
    df_intron = get_cpcdh_intron(
        cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
    ).query("name.str.lower().str.startswith('pcdha')")

    df = (
        (
            pd
            .read_feather(cfg["data_dir"] / "result" / "reads.feather")
            .assign(assemble=lambda df: df["clone"].map(clone2assemble))
            .query("assemble == @assemble")
            .reset_index(drop=True)
        )
        .assign(**{
            f"splice.{name}": lambda df, chrom=chrom, start=start, end=end: df[
                "ref_blocks"
            ].map(
                lambda ref_blocks, chrom=chrom, start=start, end=end: check.connect(
                    list(blocks_string2tuple(ref_blocks)), chrom, start, end, "+"
                )
            )
            for chrom, start, end, name in zip(
                df_intron["chrom"],
                df_intron["start"],
                df_intron["end"],
                df_intron["name"],
            )
        })
        .assign(**{
            f"precursor.{name}": lambda df, chrom=chrom, start=start: df[
                "ref_blocks"
            ].map(
                lambda ref_blocks, chrom=chrom, start=start - cfg["cover_threshold"], end=start + cfg["cover_threshold"]: (
                    check.cover(
                        list(blocks_string2tuple(ref_blocks)), chrom, start, end, "+"
                    )
                )
            )
            for chrom, start, name in zip(
                df_intron["chrom"], df_intron["start"], df_intron["name"]
            )
        })
    )

    df = (
        df
        .groupby(by=["exp", "protein", "clone", "rep", "query_name"], as_index=False)
        .agg(**{
            f"{opt}.{name}": pd.NamedAgg(column=f"{opt}.{name}", aggfunc="any")
            for name in df_intron["name"]
            for opt in ["splice", "precursor"]
        })
        .copy()
    )

    select_total_count = SelectTotalCount(cfg)
    df = (
        df
        .assign(treat=lambda df: df["clone"].map(clone2treat))
        .groupby(["exp", "protein", "treat"], as_index=False)
        .agg(**{
            f"{opt}.{name}": pd.NamedAgg(column=f"{opt}.{name}", aggfunc="sum")
            for name in df_intron["name"]
            for opt in ["splice", "precursor"]
        })
        .copy()
        .assign(
            total_count=lambda df: [
                select_total_count(exp, protein, treat)
                for exp, protein, treat in zip(df["exp"], df["protein"], df["treat"])
            ]
        )
    )

    df = (
        df
        .melt(
            id_vars=["exp", "protein", "treat", "total_count"],
            value_vars=[
                f"{opt}.{name}"
                for name in df_intron["name"]
                for opt in ["splice", "precursor"]
            ],
            var_name="opt_intron",
            value_name="count",
        )
        .assign(
            opt=lambda df: df["opt_intron"].str.split(".", expand=True)[0],
            intron=lambda df: df["opt_intron"].str.split(".", expand=True)[1],
        )
        .pivot_table(
            values="count",
            index=["exp", "protein", "treat", "total_count", "intron"],
            columns="opt",
        )
        .reset_index()
    )

    df = (
        pd
        .merge(
            df,
            df_intron[["chrom", "start", "end", "name"]],
            how="left",
            left_on="intron",
            right_on="name",
            validate="many_to_one",
        )
        .sort_values(by=["exp", "protein", "treat", "start"], ignore_index=True)
        .drop(columns="intron")
    )

    df.to_csv(cfg["data_dir"] / "result" / f"{assemble}_splice.csv", index=False)


def read_start_around_exon_start(cfg: dict, assemble: str) -> None:
    cpcdh_file = cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    tsses = df_cpcdh.query("name.str.lower().str.startswith('pcdha')")[
        ["start", "name"]
    ].reset_index(drop=True)

    df = (
        pd
        .read_feather(cfg["data_dir"] / "result" / "reads.feather")
        .assign(assemble=lambda df: df["clone"].map(clone2assemble))
        .query("assemble == @assemble")
        .assign(
            treat=lambda df: df["clone"].map(clone2treat),
            exp_protein_treat=lambda df: (
                df["exp"] + "_" + df["protein"] + "_" + df["treat"]
            ),
        )
        .assign(
            read_start=lambda df: df["ref_blocks"].map(
                lambda ref_blocks: check.start(list(blocks_string2tuple(ref_blocks)))
            )
        )[["exp_protein_treat", "read_start"]]
        .value_counts()
        .reset_index()
    )

    df = df.assign(**{
        exon: lambda df, tss=tss: df["read_start"] - tss
        for tss, exon in zip(tsses["start"], tsses["name"])
    })

    tss_extend = cfg["tss_extend"]
    df = (
        df
        .melt(
            id_vars=["exp_protein_treat", "read_start", "count"],
            value_vars=tsses["name"].tolist(),
            var_name="exon",
            value_name="relative",
        )
        .query("relative >= -@tss_extend and relative <= @tss_extend")
        .reset_index(drop=True)
        .assign(
            exon=lambda df: pd.Categorical(
                df["exon"],
                categories=tsses["name"].tolist(),
                ordered=True,
            )
        )
    )

    select_total_count = SelectTotalCount(cfg)
    pdf_files = []
    with pypdf.PdfWriter() as pdf_writer:
        for exp_protein_treat in df["exp_protein_treat"].unique():
            df_slice = (
                df
                .query("exp_protein_treat == @exp_protein_treat")
                .reset_index(drop=True)[["exon", "relative", "count"]]
                .set_index(["exon", "relative"])
                .reindex(
                    index=pd.MultiIndex.from_product(
                        [
                            tsses["name"],
                            list(range(-tss_extend, tss_extend + 1)),
                        ],
                        names=["exon", "relative"],
                    ),
                    fill_value=0,
                )
                .reset_index()
            )

            total_count = select_total_count(*exp_protein_treat.split("_"))
            df_slice = df_slice.assign(
                RPM=lambda df, total_count=total_count: (
                    df["count"] / total_count * 1_000_000
                ),
                RPM_round=lambda df: df["RPM"].round(2),
            )

            pdf_file = (
                cfg["data_dir"] / "result" / f"{assemble}_{exp_protein_treat}.pdf"
            )
            (
                ggplot(df_slice, mapping=aes(x="relative", y="exon"))
                + geom_tile(aes(fill="RPM"), color="#000000")
                + geom_text(aes(label="RPM_round"), size=6)
                + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
                + scale_y_discrete(limits=tsses["name"].tolist()[::-1])
                + theme(
                    axis_text_x=element_text(angle=90, ma="right"), figure_size=(20, 20)
                )
                + labs(title=exp_protein_treat, x="position")
            ).save(pdf_file)

            pdf_writer.append(pdf_file)
            pdf_files.append(pdf_file)

        pdf_writer.write(
            cfg["data_dir"] / "result" / f"{assemble}_read_start_around_tss.pdf"
        )

    for pdf_file in pdf_files:
        pdf_file.unlink()


def get_exon_pre(cfg: dict):
    dfs = {
        "hg19": [],
        "mm10": [],
    }
    for exp, protein, treat, bamfile in get_merge_bam(cfg).itertuples(index=False):
        bamfile = pathlib.Path(bamfile)
        assemble = treat2assemble(treat)

        df_cpcdh = (
            pd
            .read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0)
            .query("name.str.lower().str.startswith('pcdha')")
            .reset_index(drop=True)
        )
        df_cpcdh = df_cpcdh.assign(
            pre_start=lambda df: (
                [df.loc[0, "start"].item() - cfg["size_before_first"]]
                + df["end"].to_list()[:-1]
            ),
            pre_end=lambda df: df["start"],
        )

        exon_counts = []
        pre_counts = []
        with pysam.AlignmentFile(os.fspath(bamfile)) as sam:
            for chrom, start, end, pre_start, pre_end in zip(
                df_cpcdh["chrom"],
                df_cpcdh["start"],
                df_cpcdh["end"],
                df_cpcdh["pre_start"],
                df_cpcdh["pre_end"],
            ):
                exon_counts.append(
                    sum(
                        1
                        for read in sam.fetch(chrom, start, end)
                        if is_mapped_primary_first(read)
                    )
                )

                pre_counts.append(
                    sum(
                        1
                        for read in sam.fetch(chrom, pre_start, pre_end)
                        if is_mapped_primary_first(read)
                    )
                )

        df = df_cpcdh[["chrom", "start", "end", "name", "pre_start", "pre_end"]].assign(
            exon_count=exon_counts,
            pre_count=pre_counts,
            exp=exp,
            protein=protein,
            treat=treat,
        )

        dfs[assemble].append(df)

    for assemble in ["hg19", "mm10"]:
        pd.concat(dfs[assemble], ignore_index=True).to_csv(
            cfg["data_dir"] / "result" / f"{assemble}_exon_pre.csv", index=False
        )

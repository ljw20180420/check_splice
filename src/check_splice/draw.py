import os
from collections.abc import Iterable

import matplotlib
import pandas as pd
import pypdf
from plotnine import (
    aes,
    element_text,
    geom_text,
    geom_tile,
    ggplot,
    ggtitle,
    labs,
    scale_fill_gradient,
    scale_y_discrete,
    theme,
)

matplotlib.use("agg")


def splice_heatmap(cfg: dict):
    df = pd.read_csv(cfg["data_dir"] / "result" / "splice.csv", header=0)
    df = df.assign(
        name=lambda df: pd.Categorical(
            df["name"], categories=df["name"].drop_duplicates(), ordered=True
        ),
        protein_wt=lambda df: (
            df["protein"] + "_" + df["is_WT"].map({True: "control", False: "delta"})
        ),
    ).assign(**{
        "splice (RPM)": lambda df: df["connect"] / df["total_count"] * 1000_000,
        "precursor (RPM)": lambda df: df["cover.start"] / df["total_count"] * 1000_000,
        "precursor to splice": lambda df: (df["cover.start"] / df["connect"]).replace(
            float("inf"), float("nan")
        ),
    })

    targets = [
        "splice (RPM)",
        "precursor (RPM)",
        "precursor to splice",
    ]
    for target in targets:
        df = df.assign(**{
            f"{target}_round": lambda df, target=target: df[target].round(2)
        })
    with pypdf.PdfWriter() as pdf_writer:
        for exp in df["exp"].unique():
            for target in targets:
                (
                    ggplot(
                        data=df.query("name.str.startswith('PCDHA') and exp == @exp"),
                        mapping=aes(x="name", y="protein_wt"),
                    )
                    + ggtitle(exp)
                    + geom_tile(aes(fill=target), color="#000000")
                    + geom_text(aes(label=f"{target}_round"), size=6)
                    + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
                    + theme(axis_text_x=element_text(angle=90, ma="right"))
                ).save(cfg["data_dir"] / "result" / f"{target}_{exp}.pdf")

                pdf_writer.append(cfg["data_dir"] / "result" / f"{target}_{exp}.pdf")

        pdf_writer.write(cfg["data_dir"] / "result" / "splice.pdf")

        for target in targets:
            for exp in df["exp"].unique():
                (cfg["data_dir"] / "result" / f"{target}_{exp}.pdf").unlink()


def around_heatmap(
    cfg: dict,
    df: pd.DataFrame,
    center_names: Iterable[str],
    center_axis_name: str,
    slice: str,
) -> os.PathLike:
    title = f"{slice}"
    pdf_file = cfg["data_dir"] / "result" / f"{title}.pdf"

    (
        ggplot(
            df.assign(round=lambda df: df["count"].round(2)),
            mapping=aes(x="relative", y=center_axis_name),
        )
        + geom_tile(aes(fill="count"), color="#000000")
        + geom_text(aes(label="round"), size=6)
        + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
        + scale_y_discrete(limits=list(center_names)[::-1])
        + theme(axis_text_x=element_text(angle=90, ma="right"), figure_size=(20, 20))
        + labs(title=title)
    ).save(pdf_file)

    return pdf_file

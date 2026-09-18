import matplotlib
import pandas as pd
import pypdf
from plotnine import (
    aes,
    element_text,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_fill_gradient,
    theme,
)

matplotlib.use("agg")


def splice_heatmap(cfg: dict, assemble: str):
    df = pd.read_csv(cfg["data_dir"] / "result" / f"{assemble}_splice.csv", header=0)
    df = df.assign(
        name=lambda df: pd.Categorical(
            df["name"], categories=df["name"].unique(), ordered=True
        ),
        protein_treat=lambda df: df["protein"] + "_" + df["treat"],
    ).assign(**{
        "splice (RPM)": lambda df: df["splice"] / df["total_count"] * 1_000_000,
        "precursor (RPM)": lambda df: df["precursor"] / df["total_count"] * 1_000_000,
        "splice %": lambda df: df["splice"] / (df["precursor"] + df["splice"]) * 100,
    })

    targets = [
        "splice (RPM)",
        "precursor (RPM)",
        "splice %",
    ]
    for target in targets:
        df = df.assign(**{
            f"{target}_round": lambda df, target=target: df[target].round(2)
        })
    with pypdf.PdfWriter() as pdf_writer:
        for exp in df["exp"].unique():
            for target in targets:
                df_slice = df.query("exp == @exp")
                (
                    ggplot(
                        data=df_slice,
                        mapping=aes(x="name", y="protein_treat"),
                    )
                    + geom_tile(aes(fill=target), color="#000000")
                    + geom_text(aes(label=f"{target}_round"), size=6)
                    + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
                    + theme(axis_text_x=element_text(angle=90, ma="right"))
                    + labs(title=exp, x="splice", y="sample")
                ).save(cfg["data_dir"] / "result" / f"{target}_{exp}.pdf")

                pdf_writer.append(cfg["data_dir"] / "result" / f"{target}_{exp}.pdf")

        pdf_writer.write(cfg["data_dir"] / "result" / f"{assemble}_splice.pdf")

        for target in targets:
            for exp in df["exp"].unique():
                (cfg["data_dir"] / "result" / f"{target}_{exp}.pdf").unlink()

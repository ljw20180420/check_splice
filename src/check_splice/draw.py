import matplotlib
import pandas as pd
import pypdf
from plotnine import (
    aes,
    element_text,
    geom_tile,
    ggplot,
    scale_fill_gradient,
    theme,
)

matplotlib.use("agg")


def splice_heatmap(cfg):
    df = pd.read_csv(cfg["data_dir"] / "result" / "splice.csv", header=0)
    df = df.assign(
        name=lambda df: pd.Categorical(
            df["name"], categories=df["name"].drop_duplicates(), ordered=True
        ),
        exp_protein_wt=lambda df: (
            df["exp"]
            + "_"
            + df["protein"]
            + "_"
            + df["is_WT"].map({True: "wt", False: "non"})
        ),
    ).assign(**{
        "n.connect": lambda df: df["connect"] / df["total_count"],
        "n.cover.start": lambda df: df["cover.start"] / df["total_count"],
        "n.cover.end": lambda df: df["cover.end"] / df["total_count"],
        "cover.start_to_connect": lambda df: df["cover.start"] / df["connect"],
        "cover.end_to_connect": lambda df: df["cover.end"] / df["connect"],
    })

    targets = [
        "connect",
        "cover.start",
        "cover.end",
        "n.connect",
        "n.cover.start",
        "n.cover.end",
        "cover.start_to_connect",
        "cover.end_to_connect",
    ]
    with pypdf.PdfWriter() as pdf_writer:
        for target in targets:
            output = cfg["data_dir"] / "result" / f"{target}.pdf"
            (
                ggplot(data=df, mapping=aes(x="name", y="exp_protein_wt", fill=target))
                + geom_tile(color="#000000")
                + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
                + theme(axis_text_x=element_text(angle=90, ma="right"))
            ).save(output)

            pdf_writer.append(output)

        pdf_writer.write(cfg["data_dir"] / "result" / "splice.pdf")

        for target in targets:
            output = cfg["data_dir"] / "result" / f"{target}.pdf"
            output.unlink()

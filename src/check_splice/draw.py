import matplotlib
import pandas as pd
from plotnine import (
    aes,
    element_text,
    geom_raster,
    ggplot,
    scale_fill_gradient,
    theme,
)

matplotlib.use("agg")


def splice_heatmap(cfg):
    df = pd.read_csv(cfg["data_dir"] / "result" / "splice.csv", header=0)
    introns = df["name"].drop_duplicates()
    df["name"] = pd.Categorical(df["name"], categories=introns, ordered=True)
    df["exp_protein_wt"] = (
        df["exp"]
        + "_"
        + df["protein"]
        + "_"
        + df["is_WT"].map({True: "wt", False: "non"})
    )

    for target in ["connect", "cover.start", "cover.end"]:
        (
            ggplot(data=df, mapping=aes(x="name", y="exp_protein_wt", fill=target))
            + geom_raster()
            + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
            + theme(axis_text_x=element_text(angle=90, ma="right"))
        ).save(cfg["data_dir"] / "result" / f"{target}.png")

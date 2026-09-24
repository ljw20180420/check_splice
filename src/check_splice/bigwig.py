import pandas as pd

from .common import substract_bigwig, write_bigwig
from .utils import get_merge_bam, map_to_wild_type_merge, treat2assemble, treat2diff


def construct_artifact_bw(cfg: dict, assemble: str) -> None:
    df_splice = (
        pd
        .read_csv(cfg["data_dir"] / "result" / f"{assemble}_splice.csv", header=0)
        .query("name.str.lower().str.startswith('pcdha')")
        .reset_index(drop=True)
        .assign(**{
            "splice %": lambda df: df["splice"] / (df["splice"] + df["precursor"]) * 100
        })
    )

    for exp, protein, treat, bamfile in get_merge_bam(cfg).itertuples(index=False):
        df_splice_slice = df_splice.query(
            "exp == @exp and protein == @protein and treat == @treat"
        ).reset_index(drop=True)

        if len(df_splice_slice) == 0:
            continue

        (cfg["data_dir"] / "result" / "bw").mkdir(exist_ok=True, parents=True)
        bw_file = cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{treat}.bw"

        df_pv = (
            df_splice_slice[["start", "splice %"]]
            .assign(
                end=lambda df: df["start"] + 1,
                value=lambda df: df["splice %"].fillna(0.0),
            )[["start", "end", "value"]]
            .query("value > 0")
            .sort_values(by=["start"], ignore_index=True)
        )

        write_bigwig(
            chrom=cfg[assemble]["chrom"],
            chrom_size=cfg[assemble]["length"],
            starts=df_pv["start"].to_numpy(),
            ends=df_pv["end"].to_numpy(),
            values=df_pv["value"].to_numpy(),
            bigwig_file=bw_file,
        )

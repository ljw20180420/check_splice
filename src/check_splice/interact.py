import numpy as np
import pandas as pd

from .common import (
    get_cpcdh_intron,
    interact2pairs,
    pairs2bedpe,
)
from .utils import (
    SelectTotalCount,
    clone2treat,
    get_merge_bam,
)


def get_interact(cfg: dict) -> None:
    df = (
        pd
        .read_csv(cfg["data_dir"] / "result" / "expand_splice.csv", header=0)
        .rename(
            columns={
                "query_name": "name",
                "chrom1": "sourceChrom",
                "start1": "sourceStart",
                "end1": "sourceEnd",
                "strand1": "sourceStrand",
                "chrom2": "targetChrom",
                "start2": "targetStart",
                "end2": "targetEnd",
                "strand2": "targetStrand",
            }
        )
        .assign(
            chrom=lambda df: df["sourceChrom"],
            chromStart=lambda df: np.minimum(df["sourceStart"], df["targetStart"]),
            chromEnd=lambda df: np.maximum(df["sourceEnd"], df["targetEnd"]),
            match_percent=lambda df: (
                df["match_base"] / df["joint_block"].str.len() * 100
            ),
            score=lambda df: 1000 - 10 * df["match_percent"],
            value=1.0,
            exp=lambda df: (
                df["exp"] + "_" + df["protein"] + "_" + df["clone"].map(clone2treat)
            ),
            color=0,
            sourceName=".",
            targetName=".",
        )[
            [
                "chrom",
                "chromStart",
                "chromEnd",
                "name",
                "score",
                "value",
                "exp",
                "color",
                "sourceChrom",
                "sourceStart",
                "sourceEnd",
                "sourceName",
                "sourceStrand",
                "targetChrom",
                "targetStart",
                "targetEnd",
                "targetName",
                "targetStrand",
            ]
        ]
        .sort_values(by=["chrom", "chromStart", "chromEnd"], ignore_index=True)
    )

    (cfg["data_dir"] / "result" / "hic" / "interact").mkdir(exist_ok=True, parents=True)
    for exp in df["exp"].unique():
        df.query("exp == @exp").to_csv(
            cfg["data_dir"] / "result" / "hic" / "interact" / f"{exp}.bed",
            sep="\t",
            header=False,
            index=False,
        )


def interact_to_pairs(cfg: dict) -> None:
    (cfg["data_dir"] / "result" / "hic" / "pairs").mkdir(parents=True, exist_ok=True)
    df_merge = get_merge_bam(cfg)
    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        interact2pairs(
            interact_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "interact"
            / f"{exp}_{protein}_{treat}.bed",
            pairs_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "pairs"
            / f"{exp}_{protein}_{treat}.pairs",
        )


def pairs_to_bedpe(cfg: dict) -> None:
    (cfg["data_dir"] / "result" / "hic" / "bedpe").mkdir(parents=True, exist_ok=True)
    select_total_count = SelectTotalCount(cfg)
    df_merge = get_merge_bam(cfg)
    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        pairs2bedpe(
            pairs_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "pairs"
            / f"{exp}_{protein}_{treat}.pairs",
            bedpe_file=cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_{treat}.bedpe",
            total_count=select_total_count(exp, protein, treat),
        )


class BedpeJustIntronFilter:
    def __init__(self, cfg: dict) -> None:
        self.df_introns = {
            assemble: get_cpcdh_intron(
                cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
            )
            for assemble in ["hg19", "mm10"]
        }

    def __call__(
        self, starts: pd.Series, ends: pd.Series, assembles: str | pd.Series
    ) -> pd.Series:
        masks = starts <= ends
        df = pd.DataFrame({
            "start": starts.where(masks, ends),
            "end": ends.where(masks, starts),
            "assemble": assembles,
        })
        for assemble in df["assemble"].unique():
            df["just"] = (
                df
                .query("assemble == @assemble")
                .assign(**{
                    name: lambda df, start=start, end=end: (
                        (df["start"] == start) & (df["end"] == end)
                    )
                    for start, end, name in zip(
                        self.df_introns[assemble]["start"],
                        self.df_introns[assemble]["end"],
                        self.df_introns[assemble]["name"],
                    )
                })[self.df_introns[assemble]["name"].tolist()]
                .any(axis=1)
            )

        return df["just"]

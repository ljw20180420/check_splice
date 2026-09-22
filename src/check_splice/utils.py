import os

import pandas as pd


def get_sample_bam(cfg: dict) -> pd.DataFrame:
    exps = []
    proteins = []
    clones = []
    reps = []
    files = []
    for bamfile in os.listdir(cfg["data_dir"] / "bam"):
        if not bamfile.endswith(".bam"):
            continue

        exp, protein, clone, rep = bamfile.removesuffix(".bam").split("_")
        exps.append(exp)
        proteins.append(protein)
        clones.append(clone)
        reps.append(rep)
        files.append(os.fspath(cfg["data_dir"] / "bam" / bamfile))

    return pd.DataFrame({
        "exp": exps,
        "protein": proteins,
        "clone": clones,
        "rep": reps,
        "file": files,
    })


def get_merge_bam(cfg: dict) -> pd.DataFrame:
    exps = []
    proteins = []
    treats = []
    files = []
    for bamfile in os.listdir(cfg["data_dir"] / "bam" / "merge"):
        if not bamfile.endswith(".bam"):
            continue

        exp, protein, treat = bamfile.removesuffix(".bam").split("_")
        exps.append(exp)
        proteins.append(protein)
        treats.append(treat)
        files.append(os.fspath(cfg["data_dir"] / "bam" / "merge" / bamfile))

    return (
        pd
        .DataFrame({
            "exp": exps,
            "protein": proteins,
            "treat": treats,
            "file": files,
        })
        .assign(assemble=lambda df: df["treat"].map(treat2assemble))
        .sort_values(by=["assemble", "exp", "protein", "treat"], ignore_index=True)
        .drop(columns=["assemble"])
    )


def clone2assemble(clone: str) -> str:
    if clone.startswith("mm"):
        return "mm10"
    return "hg19"


def clone2treat(clone: str) -> str:
    if clone.startswith("mmWT"):
        return "mmcontrol"
    if clone.startswith("WT"):
        return "control"
    if clone.startswith("mm"):
        return "mmtreat"
    return "treat"


def treat2assemble(treat: str) -> str:
    if treat.startswith("mm"):
        return "mm10"
    return "hg19"


def treat2diff(treat: str) -> str:
    if treat.startswith("mm"):
        return "mmdiff"
    return "diff"


def blocks_string2tuple(blocks: str):
    for block in blocks.split(";"):
        chrom, start, end, strand = block.split(":")
        yield chrom, int(start), int(end), strand


class SelectTotalCount:
    def __init__(self, cfg: dict) -> None:
        self.df_total = (
            pd
            .read_csv(cfg["data_dir"] / "result" / "total_count.csv", header=0)
            .assign(treat=lambda df: df["clone"].map(clone2treat))
            .groupby(by=["exp", "protein", "treat"], as_index=False)["total_count"]
            .sum()
        )

    def __call__(self, exp: str, protein: str, treat: str) -> int:
        return self.df_total.query(
            "exp == @exp and protein == @protein and treat == @treat"
        )["total_count"].item()


def map_to_wild_type_merge(exp: str, protein: str, treat: str) -> tuple[str, str, str]:
    if protein == "merge":
        wt_protein = "merge"
    elif exp == "clip":
        wt_protein = "WT"
    else:
        wt_protein = protein
    wt_treat = "mmcontrol" if treat.startswith("mm") else "control"

    return exp, wt_protein, wt_treat

import os

import matplotlib.pyplot as plt
import pandas as pd
import pypdf
from coolbox.api import *


def get_total_count(cfg: dict, exp_protein_wt: str) -> int:
    df_total = pd.read_csv(cfg["data_dir"] / "result" / "total_count.csv", header=0)
    df_total = (
        df_total
        .assign(
            exp_protein_wt=lambda df: (
                df["exp"]
                + "_"
                + df["protein"]
                + "_"
                + df["clone"].map(
                    lambda ele: "control" if ele.startswith("WT") else "delta"
                )
            )
        )
        .groupby("exp_protein_wt")["total_count"]
        .sum()
        .reset_index()
    )

    if exp_protein_wt:
        total_count = df_total.loc[
            df_total["exp_protein_wt"] == exp_protein_wt, "total_count"
        ].item()
    else:
        total_count = df_total["total_count"].sum()

    return total_count


def pairs_to_bedpe(cfg: dict) -> None:
    (cfg["data_dir"] / "result" / "hic" / "bedpe").mkdir(parents=True, exist_ok=True)
    for pairs_file in os.listdir(cfg["data_dir"] / "result" / "hic" / "pairs"):
        if pairs_file == "ff.pairs" or pairs_file == "rr.pairs":
            exp_protein_wt = ""
        else:
            exp_protein_wt = pairs_file.rsplit("_", 1)[0]
        total_count = get_total_count(cfg, exp_protein_wt)

        pairs_file = cfg["data_dir"] / "result" / "hic" / "pairs" / pairs_file
        bedpe_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / pairs_file.with_suffix(".bedpe").name
        )

        df = pd.read_csv(
            pairs_file,
            sep="\t",
            skiprows=1,
            names=[
                "readID",
                "chrom1",
                "pos1",
                "chrom2",
                "pos2",
                "strand1",
                "strand2",
            ],
        )

        df = (
            df
            .rename(
                columns={
                    "pos1": "end1",
                    "pos2": "end2",
                }
            )
            .assign(
                start1=lambda df: df["end1"] - 1,
                start2=lambda df: df["end2"] - 1,
            )
            .groupby([
                "chrom1",
                "start1",
                "end1",
                "strand1",
                "chrom2",
                "start2",
                "end2",
                "strand2",
            ])
            .agg(
                name=pd.NamedAgg("readID", "first"),
                score=pd.NamedAgg("readID", "count"),
            )
            .reset_index()
            .assign(
                score=lambda df, total_count=total_count: (
                    df["score"] / total_count * 1000000
                )
            )[
                [
                    "chrom1",
                    "start1",
                    "end1",
                    "chrom2",
                    "start2",
                    "end2",
                    "name",
                    "score",
                    "strand1",
                    "strand2",
                ]
            ]
        )

        df.to_csv(bedpe_file, sep="\t", header=False, index=False)


def draw_links(
    cfg: dict,
    exp_protein_wt: str,
    cluster: str,
):
    ff_bedpe_file = (
        cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp_protein_wt}_ff.bedpe"
    )
    rr_bedpe_file = (
        cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp_protein_wt}_rr.bedpe"
    )

    protein_wt = exp_protein_wt.split("_", 1)[1]
    protein = protein_wt.rsplit("_")[0]
    color = cfg["color"][protein]
    frame = (
        XAxis(name="hg19")
        + BEDPE(os.fspath(ff_bedpe_file))
        + Color(color)
        + TrackHeight(5)
        + Title("splice")
        + BED(
            os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
            display="collapsed",
            labels=False,
        )
        + Title("gene")
        + BEDPE(os.fspath(rr_bedpe_file))
        + Color(color)
        + Inverted()
        + TrackHeight(5)
        + Title("splice")
        + FrameTitle(f"{protein_wt}:{cluster}")
    )
    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)
    link_file = (
        cfg["data_dir"] / "result" / "hic" / "draw" / f"{exp_protein_wt}_{cluster}.pdf"
    )
    region = cfg[cluster]
    fig = frame.plot(f"{region['chrom']}:{region['start']}-{region['end']}")
    fig.savefig(os.fspath(link_file))
    plt.close(fig)

    return link_file


def get_sorted_exp_protein_cluster_wts(cfg: dict) -> list[str]:
    exp_protein_wts = []
    for bedpe_file in os.listdir(cfg["data_dir"] / "result" / "hic" / "bedpe"):
        if not bedpe_file.endswith("ff.bedpe"):
            continue
        exp_protein_wt = bedpe_file.removesuffix("ff.bedpe").rstrip("_")
        if not exp_protein_wt:
            continue
        exp_protein_wts.append(exp_protein_wt)

    return (
        pd
        .DataFrame({"exp_protein_wt": exp_protein_wts})["exp_protein_wt"]
        .str.split("_", expand=True)
        .rename(columns={0: "exp", 1: "protein", 2: "wt"})
        .assign(cluster=lambda df: [["alpha", "beta", "gamma"]] * len(df))
        .explode(column="cluster", ignore_index=True)[
            ["exp", "protein", "cluster", "wt"]
        ]
        .sort_values(by=["exp", "protein", "cluster", "wt"], ignore_index=True)
    )


def draw_links_exp(cfg: dict, exp: str):
    pdf_files = []
    df_sort = get_sorted_exp_protein_cluster_wts(cfg)
    with pypdf.PdfWriter() as pdf_writer:
        for experi, protein, cluster, wt in df_sort.itertuples(index=False):
            if experi != exp:
                continue
            exp_protein_wt = f"{experi}_{protein}_{wt}"
            pdf_file = draw_links(
                cfg,
                exp_protein_wt,
                cluster,
            )
            pdf_writer.append(pdf_file)
            pdf_files.append(pdf_file)

        pdf_writer.write(
            cfg["data_dir"] / "result" / "hic" / "draw" / f"{exp}_links.pdf"
        )

    for pdf_file in pdf_files:
        pdf_file.unlink()


def draw_links_all(cfg: dict):
    draw_links_exp(cfg, "total")
    draw_links_exp(cfg, "rna")
    draw_links_exp(cfg, "pro")
    draw_links_exp(cfg, "clip")

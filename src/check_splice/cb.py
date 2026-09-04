import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import oxbow as ox
import pandas as pd
import pyBigWig
import pypdf
from coolbox.api import *
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer.compute_features_levels import compute_features_levels

from .sam import get_precursor_pos


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
    shutil.rmtree(cfg["data_dir"] / "result" / "hic" / "bedpe", ignore_errors=True)
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
                    df["score"] / total_count * 1_000_000
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


def diff_bedpe(cfg: dict, exp: str, protein: str, orientation: str) -> None:
    treat_file = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_delta_{orientation}.bedpe"
    )
    if exp != "clip":
        control_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_control_{orientation}.bedpe"
        )
    else:
        control_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_WT_control_{orientation}.bedpe"
        )

    df_treat = pd.read_csv(
        treat_file,
        sep="\t",
        names=[
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
        ],
    )
    df_control = pd.read_csv(
        control_file,
        sep="\t",
        names=[
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
        ],
    )
    df = df_treat.merge(
        df_control,
        on=["chrom1", "start1", "end1", "chrom2", "start2", "end2"],
        how="outer",
    ).assign(
        name=lambda df: df["name_x"].combine_first(df["name_y"]),
        score=lambda df: df["score_x"].fillna(0) - df["score_y"].fillna(0),
        strand1=lambda df: df["strand1_x"].combine_first(df["strand1_y"]),
        strand2=lambda df: df["strand2_x"].combine_first(df["strand2_y"]),
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

    df.query("score > 0").to_csv(
        (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_diff_{orientation}_increase.bedpe"
        ),
        sep="\t",
        index=False,
        header=False,
    )
    df.query("score < 0").assign(score=lambda df: -df["score"]).to_csv(
        (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_diff_{orientation}_decrease.bedpe"
        ),
        sep="\t",
        index=False,
        header=False,
    )


def diff_bedpe_all(cfg: dict) -> None:
    pairs_to_bedpe(cfg)
    for exp in ["total", "rna", "pro", "clip"]:
        for protein in ["NP220", "MPP8", "PPHLN1", "TASOR"]:
            for orientation in ["ff", "rr"]:
                diff_bedpe(cfg, exp, protein, orientation)


def draw_links(
    cfg: dict,
    exp: str,
    protein: str,
    cluster: str,
) -> os.PathLike:
    if exp != "clip":
        control_f = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_control_ff.bedpe"
        )
        control_r = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_control_rr.bedpe"
        )
    else:
        control_f = (
            cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp}_WT_control_ff.bedpe"
        )
        control_r = (
            cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp}_WT_control_rr.bedpe"
        )

    treat_f = (
        cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp}_{protein}_delta_ff.bedpe"
    )
    treat_r = (
        cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp}_{protein}_delta_rr.bedpe"
    )

    diff_f_increase = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_diff_ff_increase.bedpe"
    )
    diff_r_increase = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_diff_rr_increase.bedpe"
    )
    diff_f_decrease = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_diff_ff_decrease.bedpe"
    )
    diff_r_decrease = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_diff_rr_decrease.bedpe"
    )

    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)
    treat = "delta" if exp != "clip" else "tag"

    frame = (
        XAxis(name="hg19")
        + BEDPE(
            os.fspath(control_f),
            color=cfg["color"]["WT"],
            height=5,
            title="control",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BEDPE(
            os.fspath(control_r),
            color=cfg["color"]["WT"],
            height=5,
            title="control",
            orientation="inverted",
        )
        + BEDPE(
            os.fspath(treat_f),
            color=cfg["color"][protein],
            height=5,
            title=treat,
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BEDPE(
            os.fspath(treat_r),
            color=cfg["color"][protein],
            height=5,
            title=treat,
            orientation="inverted",
        )
        + BEDPE(
            os.fspath(diff_f_increase),
            color=cfg["color"]["INCREASE"],
            height=5,
            title="diff",
        )
        + BEDPECoverage(os.fspath(diff_f_decrease), color=cfg["color"]["DECREASE"])
        + BED(
            os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BEDPE(
            os.fspath(diff_r_increase),
            color=cfg["color"]["INCREASE"],
            height=5,
            title="diff",
        )
        + BEDPECoverage(
            os.fspath(diff_r_decrease),
            color=cfg["color"]["DECREASE"],
            orientation="inverted",
        )
        + FrameTitle(protein)
    )
    link_file = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "draw"
        / f"{exp}_{protein}_{cluster}_links.pdf"
    )
    chrom = cfg[cluster]["chrom"]
    start = cfg[cluster]["start"]
    end = cfg[cluster]["end"]
    fig = frame.plot(chrom, start, end)
    fig.savefig(os.fspath(link_file))
    plt.close(fig)

    return link_file


def draw_covers(
    cfg: dict,
    exp: str,
    protein: str,
    cluster: str,
) -> os.PathLike:
    if exp != "clip":
        control_f = (
            cfg["data_dir"]
            / "bam"
            / "precursor"
            / "merge"
            / f"{exp}_{protein}_control.f.bw"
        )
        control_r = (
            cfg["data_dir"]
            / "bam"
            / "precursor"
            / "merge"
            / f"{exp}_{protein}_control.r.bw"
        )
    else:
        control_f = (
            cfg["data_dir"] / "bam" / "precursor" / "merge" / f"{exp}_WT_control.f.bw"
        )
        control_r = (
            cfg["data_dir"] / "bam" / "precursor" / "merge" / f"{exp}_WT_control.r.bw"
        )

    treat = "delta" if exp != "clip" else "tag"

    treat_f = (
        cfg["data_dir"]
        / "bam"
        / "precursor"
        / "merge"
        / f"{exp}_{protein}_{treat}.f.bw"
    )
    treat_r = (
        cfg["data_dir"]
        / "bam"
        / "precursor"
        / "merge"
        / f"{exp}_{protein}_{treat}.r.bw"
    )

    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)

    chrom = cfg[cluster]["chrom"]
    start = cfg[cluster]["start"]
    end = cfg[cluster]["end"]

    max_heights = []
    for bwfile in [control_f, control_r, treat_f, treat_r]:
        with pyBigWig.open(os.fspath(bwfile)) as bw:
            max_height = bw.stats(chrom, start, end, type="max")[0]
            if max_height is not None:
                max_heights.append(max_height)

    if max_heights:
        yup = max(max_heights) * 1.1
    else:
        yup = 0

    frame = (
        XAxis(name="hg19")
        + BigWig(
            os.fspath(control_f),
            min_value=0,
            max_value=yup,
            color=cfg["color"]["WT"],
            height=5,
            title="control",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BigWig(
            os.fspath(control_r),
            min_value=0,
            max_value=yup,
            color=cfg["color"]["WT"],
            height=5,
            title="control",
            orientation="inverted",
        )
        + BigWig(
            os.fspath(treat_f),
            min_value=0,
            max_value=yup,
            color=cfg["color"][protein],
            height=5,
            title=treat,
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BigWig(
            os.fspath(treat_r),
            min_value=0,
            max_value=yup,
            color=cfg["color"][protein],
            height=5,
            title=treat,
            orientation="inverted",
        )
        + FrameTitle(protein)
    )
    link_file = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "draw"
        / f"{exp}_{protein}_{cluster}_covers.pdf"
    )
    fig = frame.plot(chrom, start, end)
    fig.savefig(os.fspath(link_file))
    plt.close(fig)

    return link_file


def draw_all(cfg: dict):
    for exp in ["total", "rna", "pro", "clip"]:
        pdf_files = []
        with pypdf.PdfWriter() as pdf_writer:
            for protein in ["NP220", "MPP8", "PPHLN1", "TASOR"]:
                for cluster in ["alpha", "beta", "gamma"]:
                    pdf_file = draw_links(
                        cfg,
                        exp=exp,
                        protein=protein,
                        cluster=cluster,
                    )
                    pdf_writer.append(pdf_file)
                    pdf_files.append(pdf_file)

                    pdf_file = draw_covers(
                        cfg,
                        exp=exp,
                        protein=protein,
                        cluster=cluster,
                    )
                    pdf_writer.append(pdf_file)
                    pdf_files.append(pdf_file)

            pdf_writer.write(cfg["data_dir"] / "result" / "hic" / "draw" / f"{exp}.pdf")

        for pdf_file in pdf_files:
            pdf_file.unlink()


def estimate_height(bam_file: os.PathLike, chrom: str, start: int, end: int) -> float:
    ds = ox.from_bam(bam_file)
    sub = ds.regions(f"{chrom}:{start}-{end}")
    df = sub.pd()
    rev_flag = np.bitwise_and(df["flag"], 0b10000) != 0
    features = []
    for idx, row in df.iterrows():
        start = row["pos"] - start
        end = row["pos"] + len(row["seq"]) - start
        strand = -1 if rev_flag.iloc[idx] else 1
        gf = GraphicFeature(
            start=start,
            end=end,
            strand=strand,
        )
        features.append(gf)
    record = GraphicRecord(sequence_length=end - start, features=features)

    feature_levels = compute_features_levels(record.features)
    if feature_levels:
        max_track_level = max(feature_levels.values())
        total_tracks = max_track_level + 1  # 0-indexed levels
    else:
        total_tracks = 1  # Fallback if no features
    calculated_height = 1.5 + (total_tracks * 0.4)

    return calculated_height


def draw_reads(
    cfg: dict,
    exp: str,
    protein: str,
):
    if exp != "clip":
        control_f = (
            cfg["data_dir"]
            / "bam"
            / "precursor"
            / "merge"
            / f"{exp}_{protein}_control.f.unify.bam"
        )
        control_r = (
            cfg["data_dir"]
            / "bam"
            / "precursor"
            / "merge"
            / f"{exp}_{protein}_control.r.unify.bam"
        )
    else:
        control_f = (
            cfg["data_dir"]
            / "bam"
            / "precursor"
            / "merge"
            / f"{exp}_WT_control.f.unify.bam"
        )
        control_r = (
            cfg["data_dir"]
            / "bam"
            / "precursor"
            / "merge"
            / f"{exp}_WT_control.r.unify.bam"
        )

    treat = "delta" if exp != "clip" else "tag"

    treat_f = (
        cfg["data_dir"]
        / "bam"
        / "precursor"
        / "merge"
        / f"{exp}_{protein}_{treat}.f.unify.bam"
    )
    treat_r = (
        cfg["data_dir"]
        / "bam"
        / "precursor"
        / "merge"
        / f"{exp}_{protein}_{treat}.r.unify.bam"
    )

    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)
    df_se = get_precursor_pos(cfg)
    for name, se, pos in zip(df_se["name"], df_se["se"], df_se["pos"]):
        chrom = cfg["chrom"]
        start = pos - 150
        end = pos + 150
        frame = (
            XAxis(name="hg19")
            + BAM(
                os.fspath(control_f),
                length_ratio_thresh=1e-5,
                color=cfg["color"]["WT"],
                height=estimate_height(control_f, chrom, start, end),
                title="control",
            )
            + BED(
                os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
                display="collapsed",
                labels=False,
                title=f"{name}:{se}",
            )
            + BAM(
                os.fspath(control_r),
                length_ratio_thresh=1e-5,
                color=cfg["color"]["WT"],
                height=estimate_height(control_r, chrom, start, end),
                title="control",
                orientation="inverted",
            )
            + BAM(
                os.fspath(treat_f),
                length_ratio_thresh=1e-5,
                color=cfg["color"][protein],
                height=estimate_height(treat_f, chrom, start, end),
                title=treat,
            )
            + BED(
                os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
                display="collapsed",
                labels=False,
                title=f"{name}:{se}",
            )
            + BAM(
                os.fspath(treat_r),
                length_ratio_thresh=1e-5,
                color=cfg["color"][protein],
                height=estimate_height(treat_r, chrom, start, end),
                title=treat,
                orientation="inverted",
            )
            + FrameTitle(protein)
        )
        link_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "draw"
            / f"{exp}_{protein}_{name}_{se}_reads.pdf"
        )
        fig = frame.plot(chrom, start, end)
        fig.savefig(os.fspath(link_file))
        plt.close(fig)

        yield link_file


def draw_reads_all(cfg: dict):
    for exp in ["total", "rna", "pro", "clip"]:
        pdf_files = []
        with pypdf.PdfWriter() as pdf_writer:
            for protein in ["NP220", "MPP8", "PPHLN1", "TASOR"]:
                for pdf_file in draw_reads(cfg, exp, protein):
                    pdf_writer.append(pdf_file)
                    pdf_files.append(pdf_file)

            pdf_writer.write(
                cfg["data_dir"] / "result" / "hic" / "draw" / f"{exp}_reads.pdf"
            )

        for pdf_file in pdf_files:
            pdf_file.unlink()

import os

import matplotlib.pyplot as plt
import numpy as np
import oxbow as ox
import pyBigWig
from coolbox.api import *
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer.compute_features_levels import compute_features_levels

from .common import (
    get_cpcdh_intron,
    read_bedpe,
)
from .utils import get_merge_bam, map_to_wild_type_merge, treat2assemble, treat2diff


def draw_link(
    cfg: dict,
    assemble: str,
    cluster: str,
    title: str,
    control_f: os.PathLike,
    treat_f: os.PathLike,
    out_f: os.PathLike,
):
    chrom = cfg[assemble][cluster]["chrom"]
    start = cfg[assemble][cluster]["start"]
    end = cfg[assemble][cluster]["end"]

    max_score = 0
    for bedpe_file in [control_f, treat_f]:
        df_bedpe = read_bedpe(bedpe_file)
        max_score = max(max_score, float(df_bedpe["score"].max()))

    diameter_to_height = f"0.5 * max_height * diameter / ({end} - {start})"
    height = 1.5
    cmap = "gray"
    frame = (
        Frame(width=18)
        + XAxis(name=assemble)
        + BEDPE(
            os.fspath(control_f),
            line_width=1,
            cmap=cmap,
            vmin=0,
            vmax=max_score,
            diameter_to_height=diameter_to_height,
            height=height,
            title="control",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BEDPE(
            os.fspath(treat_f),
            line_width=1,
            cmap=cmap,
            vmin=0,
            vmax=max_score,
            diameter_to_height=diameter_to_height,
            height=height,
            title="treat",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + FrameTitle(title)
    )

    fig = frame.plot(chrom, start, end)
    fig.savefig(os.fspath(out_f))
    plt.close(fig)


def draw_links(
    cfg: dict,
    cluster: str,
):
    df_merge = (
        get_merge_bam(cfg).query("treat.str.endswith('treat')").reset_index(drop=True)
    )
    for exp, treat in (
        df_merge[["exp", "treat"]].drop_duplicates().itertuples(index=False)
    ):
        assemble = treat2assemble(treat)
        df_merge_slice = df_merge.query("exp == @exp and treat == @treat").reset_index(
            drop=True
        )
        for protein in ["merge"] + df_merge_slice["protein"].tolist():
            treat_f = (
                cfg["data_dir"]
                / "result"
                / "hic"
                / "bedpe"
                / "cpcdh"
                / f"{exp}_{protein}_{treat}.f.bedpe"
            )

            _, wt_protein, control = map_to_wild_type_merge(exp, protein, treat)
            control_f = (
                cfg["data_dir"]
                / "result"
                / "hic"
                / "bedpe"
                / "cpcdh"
                / f"{exp}_{wt_protein}_{control}.f.bedpe"
            )

            (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(
                parents=True, exist_ok=True
            )

            link_file = (
                cfg["data_dir"]
                / "result"
                / "hic"
                / "draw"
                / f"{assemble}_{exp}_{protein}_{cluster}_links.pdf"
            )
            draw_link(
                cfg=cfg,
                assemble=assemble,
                cluster=cluster,
                title=f"{assemble}_{exp}_{protein}_{cluster}",
                control_f=control_f,
                treat_f=treat_f,
                out_f=link_file,
            )

            yield link_file


def draw_pre_exon(
    cfg: dict,
    assemble: str,
    cluster: str,
    protein: str,
    title: str,
    control_f: os.PathLike,
    treat_f: os.PathLike,
    diff_up_f: os.PathLike,
    diff_down_f: os.PathLike,
    out_f: os.PathLike,
) -> os.PathLike:
    chrom = cfg[assemble][cluster]["chrom"]
    start = cfg[assemble][cluster]["start"]
    end = cfg[assemble][cluster]["end"]

    max_heights = []
    for bwfile in [control_f, treat_f]:
        with pyBigWig.open(os.fspath(bwfile)) as bw:
            max_height = bw.stats(chrom, start, end, type="max")[0]
            if max_height is not None:
                max_heights.append(max_height)

    if max_heights:
        yup = min(max(max_heights) * 1.1, 100)
    else:
        yup = 0

    height = 1
    frame = (
        Frame(width=18)
        + XAxis(name=assemble)
        + BigWig(
            os.fspath(control_f),
            min_value=0,
            max_value=yup,
            threshold=0,
            threshold_color=cfg["color"]["WT"],
            height=height,
            title="control",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BigWig(
            os.fspath(treat_f),
            min_value=0,
            max_value=yup,
            threshold=0,
            threshold_color=cfg["color"][protein],
            height=height,
            title="treat",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BigWig(
            os.fspath(diff_up_f),
            min_value=0,
            max_value=yup,
            threshold=0,
            threshold_color=cfg["color"]["INCREASE"],
            height=height,
            title="increase",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BigWig(
            os.fspath(diff_down_f),
            min_value=0,
            max_value=yup,
            threshold=0,
            threshold_color=cfg["color"]["DECREASE"],
            height=height,
            title="decrease",
            orientation="inverted",
        )
        + FrameTitle(title)
    )
    fig = frame.plot(chrom, start, end)
    fig.savefig(os.fspath(out_f))
    plt.close(fig)


def draw_pre_exons(cfg: dict, cluster: str):
    df_merge = (
        get_merge_bam(cfg).query("treat.str.endswith('treat')").reset_index(drop=True)
    )
    for exp, protein, treat in (
        df_merge[["exp", "protein", "treat"]].drop_duplicates().itertuples(index=False)
    ):
        assemble = treat2assemble(treat)

        treat_f = cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{treat}.bw"
        _, wt_protein, control = map_to_wild_type_merge(exp, protein, treat)
        control_f = (
            cfg["data_dir"] / "result" / "bw" / f"{exp}_{wt_protein}_{control}.bw"
        )
        diff_up_f = (
            cfg["data_dir"]
            / "result"
            / "bw"
            / f"{exp}_{protein}_{treat2diff(treat)}.up.bw"
        )
        diff_down_f = (
            cfg["data_dir"]
            / "result"
            / "bw"
            / f"{exp}_{protein}_{treat2diff(treat)}.down.bw"
        )

        (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)

        pre_exon_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "draw"
            / f"{assemble}_{exp}_{protein}_{cluster}_pre_exon.pdf"
        )

        draw_pre_exon(
            cfg=cfg,
            assemble=assemble,
            cluster=cluster,
            protein=protein,
            title=f"{assemble}_{exp}_{protein}_{cluster}",
            control_f=control_f,
            treat_f=treat_f,
            diff_up_f=diff_up_f,
            diff_down_f=diff_down_f,
            out_f=pre_exon_file,
        )

        yield pre_exon_file


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


def draw_read(
    cfg: dict,
    assemble: str,
    pos: int,
    protein: str,
    exon: str,
    se: str,
    title: str,
    control_f: os.PathLike,
    treat_f: os.PathLike,
    out_f: os.PathLike,
):
    chrom = cfg[assemble]["chrom"]
    start = pos - cfg["read_length"]
    end = pos + cfg["read_length"]
    frame = (
        XAxis(name=assemble)
        + BAM(
            os.fspath(control_f.with_suffix(".f.bam")),
            length_ratio_thresh=0,
            color=cfg["color"]["WT"],
            height=estimate_height(control_f.with_suffix(".f.bam"), chrom, start, end),
            title="control",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=f"{exon}:{se}",
        )
        + BAM(
            os.fspath(control_f.with_suffix(".r.bam")),
            length_ratio_thresh=0,
            color=cfg["color"]["WT"],
            height=estimate_height(control_f.with_suffix(".r.bam"), chrom, start, end),
            title="control",
            orientation="inverted",
        )
        + BAM(
            os.fspath(treat_f.with_suffix(".f.bam")),
            length_ratio_thresh=0,
            color=cfg["color"][protein],
            height=estimate_height(treat_f.with_suffix(".f.bam"), chrom, start, end),
            title="treat",
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=f"{exon}:{se}",
        )
        + BAM(
            os.fspath(treat_f.with_suffix(".r.bam")),
            length_ratio_thresh=0,
            color=cfg["color"][protein],
            height=estimate_height(treat_f.with_suffix(".r.bam"), chrom, start, end),
            title="treat",
            orientation="inverted",
        )
        + FrameTitle(title)
    )
    fig = frame.plot(chrom, start, end)
    fig.savefig(os.fspath(out_f))
    plt.close(fig)


def draw_reads(cfg: dict):
    df_merge = (
        get_merge_bam(cfg).query("treat.str.endswith('treat')").reset_index(drop=True)
    )
    for exp, protein, treat in (
        df_merge[["exp", "protein", "treat"]].drop_duplicates().itertuples(index=False)
    ):
        assemble = treat2assemble(treat)

        treat_f = (
            cfg["data_dir"]
            / "bam"
            / "merge"
            / "precursor"
            / f"{exp}_{protein}_{treat}.bam"
        )
        _, wt_protein, control = map_to_wild_type_merge(exp, protein, treat)
        control_f = (
            cfg["data_dir"]
            / "bam"
            / "merge"
            / "precursor"
            / f"{exp}_{wt_protein}_{control}.bam"
        )

        (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)
        df_se = (
            get_cpcdh_intron(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv")
            .melt(
                id_vars=["chrom", "name"],
                value_vars=["start", "end"],
                var_name="se",
                value_name="pos",
            )
            .query("name.str.lower().str.startswith('pcdha')")
            .reset_index(drop=True)
            .sort_values(by="pos", ignore_index=True)
        )

        for exon, se, pos in zip(df_se["name"], df_se["se"], df_se["pos"]):
            read_file = (
                cfg["data_dir"]
                / "result"
                / "hic"
                / "draw"
                / f"{exp}_{protein}_{exon}_{se}_reads.pdf"
            )

            draw_read(
                cfg=cfg,
                assemble=assemble,
                pos=pos,
                protein=protein,
                exon=exon,
                se=se,
                title=f"{assemble}_{exp}_{protein}_{exon}_{se}",
                control_f=control_f,
                treat_f=treat_f,
                out_f=read_file,
            )

            yield read_file

import os
import pathlib
import tempfile

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import oxbow as ox
import pandas as pd
import pyBigWig
from coolbox.api import *
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer.compute_features_levels import compute_features_levels

from .common import (
    bedpe_in_range_with_strand,
    draw_color_bar,
    interact2bedpe,
    read_interact,
    substract_bigwig,
    summation_bedpe,
)
from .interact import BedpeJustIntronFilter
from .utils import (
    SelectTotalCount,
    get_merge_bam,
    map_to_wild_type_merge,
    treat2assemble,
    treat2diff,
)


def draw_link(
    cfg: dict,
    assemble: str,
    cluster: str,
    title: str,
    df_control: pd.DataFrame,
    df_treat: pd.DataFrame,
    out_f: os.PathLike,
    bar_f: os.PathLike,
):
    chrom = cfg[assemble][cluster]["chrom"]
    start = cfg[assemble][cluster]["start"]
    end = cfg[assemble][cluster]["end"]
    df_control = bedpe_in_range_with_strand(df_control, chrom, start, end)
    df_treat = bedpe_in_range_with_strand(df_treat, chrom, start, end)

    max_score = 0
    for df_bedpe in [df_control, df_treat]:
        max_score = max(max_score, float(df_bedpe["score"].max()))

    diameter_to_height = f"0.5 * max_height * diameter / ({end} - {start})"
    height = 1.5

    cmap = "truncated_gray_r"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        control_f = tmpdir / "control.f.bedpe"
        df_control.query("strand1 == '+' and strand2 == '+'").to_csv(
            control_f, sep="\t", index=False, header=False
        )
        control_r = tmpdir / "control.r.bedpe"
        df_control.query("strand1 == '-' and strand2 == '-'").to_csv(
            control_r, sep="\t", index=False, header=False
        )
        treat_f = tmpdir / "treat_f.bedpe"
        df_treat.query("strand1 == '+' and strand2 == '+'").to_csv(
            treat_f, sep="\t", index=False, header=False
        )
        treat_r = tmpdir / "treat_r.bedpe"
        df_treat.query("strand1 == '-' and strand2 == '-'").to_csv(
            treat_r, sep="\t", index=False, header=False
        )

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
                os.fspath(control_r),
                line_width=1,
                cmap=cmap,
                vmin=0,
                vmax=max_score,
                diameter_to_height=diameter_to_height,
                height=height,
                title="control",
                orientation="inverted",
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
            + BEDPE(
                os.fspath(treat_r),
                line_width=1,
                cmap=cmap,
                vmin=0,
                vmax=max_score,
                diameter_to_height=diameter_to_height,
                height=height,
                title="treat",
                orientation="inverted",
            )
            + FrameTitle(title)
        )

        fig = frame.plot(chrom, start, end)
        fig.savefig(os.fspath(out_f))

        draw_color_bar(
            cmap=cmap,
            vmin=0,
            vmax=max_score,
            label="RPM",
            outfile=os.fspath(bar_f),
        )


def process_interact(
    cfg: dict, exp: str, protein: str, treat: str, blat: bool, just: bool
):
    _, wt_protein, control = map_to_wild_type_merge(exp, protein, treat)
    for protein_, treat_ in [(protein, treat), (wt_protein, control)]:
        df_interact = read_interact(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "interact"
            / f"{exp}_{protein_}_{treat_}.bed"
        )
        if blat:
            max_match_percent = cfg["max_match_percent"]
            df_interact = df_interact.query(
                "(1000 - score) / 10 <= @max_match_percent"
            ).reset_index(drop=True)

        select_total_count = SelectTotalCount(cfg)
        total_count = select_total_count(exp, protein_, treat_)
        df_bedpe = interact2bedpe(df_interact, total_count)
        if just:
            bedpe_just_intron_filter = BedpeJustIntronFilter(cfg)
            df_bedpe = bedpe_just_intron_filter(
                df_bedpe.query("strand1 == '+' and strand2 == '+'").reset_index(
                    drop=True
                ),
                assemble=treat2assemble(treat),
            )

        yield df_bedpe, total_count


def draw_links(cfg: dict, cluster: str, blat: bool, just: bool):
    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)
    plt.colormaps.register(
        name="truncated_gray_r",
        cmap=mcolors.LinearSegmentedColormap.from_list(
            "truncated_gray_r", plt.colormaps["gray_r"](np.linspace(0.15, 1.0, 256))
        ),
        force=True,
    )
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

        df_treats = []
        df_controls = []
        treat_total_counts = []
        control_total_counts = []
        for protein in df_merge_slice["protein"]:
            (df_treat, treat_total_count), (df_control, control_total_count) = list(
                process_interact(cfg, exp, protein, treat, blat, just)
            )

            title = f"{exp}_{protein}_{treat}_{cluster}"
            out_f = cfg["data_dir"] / "result" / "hic" / "draw" / f"{title}.links.pdf"
            bar_f = out_f.with_suffix(".colorbar.pdf")
            draw_link(
                cfg=cfg,
                assemble=assemble,
                cluster=cluster,
                title=title,
                df_control=df_control,
                df_treat=df_treat,
                out_f=out_f,
                bar_f=bar_f,
            )
            yield out_f
            yield bar_f

            df_treats.append(df_treat)
            df_controls.append(df_control)
            treat_total_counts.append(treat_total_count)
            control_total_counts.append(control_total_count)

        title = f"{exp}_summation_{treat}_{cluster}"
        out_f = cfg["data_dir"] / "result" / "hic" / "draw" / f"{title}.links.pdf"
        bar_f = out_f.with_suffix(".colorbar.links.pdf")
        draw_link(
            cfg=cfg,
            assemble=assemble,
            cluster=cluster,
            title=title,
            df_control=summation_bedpe(
                df_bedpes=df_controls, total_counts=control_total_counts
            ),
            df_treat=summation_bedpe(
                df_bedpes=df_treats, total_counts=treat_total_counts
            ),
            out_f=out_f,
            bar_f=bar_f,
        )

        yield out_f
        yield bar_f


def draw_pre_exon(
    cfg: dict,
    assemble: str,
    cluster: str,
    protein: str,
    title: str,
    control_f: os.PathLike,
    treat_f: os.PathLike,
    out_f: os.PathLike,
) -> os.PathLike:
    chrom = cfg[assemble][cluster]["chrom"]
    chrom_size = cfg[assemble]["length"]
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
    with tempfile.TemporaryDirectory() as tmpdir:
        diff_f = pathlib.Path(tmpdir) / "diff.bw"
        substract_bigwig(
            bigwig_file1=treat_f,
            bigwig_file2=control_f,
            bigwig_file_diff=diff_f,
            chrom=chrom,
            chrom_size=chrom_size,
            start=start,
            end=end,
        )

        frame = (
            Frame(width=18, margins={"left": 0.1, "right": 0.92, "bottom": 0, "top": 1})
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
                os.fspath(diff_f),
                min_value=-yup,
                max_value=yup,
                threshold=0,
                threshold_color=cfg["color"]["INCREASE"],
                color=cfg["color"]["DECREASE"],
                height=height,
                title="diff",
                spine=0.5,
            )
            + FrameTitle(title)
        )
        fig = frame.plot(chrom, start, end)
        fig.savefig(os.fspath(out_f))


def draw_pre_exons(cfg: dict, cluster: str):
    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)
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

        title = f"{assemble}_{exp}_{protein}_{cluster}"
        out_f = cfg["data_dir"] / "result" / "hic" / "draw" / f"{title}.pre_exon.pdf"
        draw_pre_exon(
            cfg=cfg,
            assemble=assemble,
            cluster=cluster,
            protein=protein,
            title=title,
            control_f=control_f,
            treat_f=treat_f,
            out_f=out_f,
        )

        yield out_f


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
            pd
            .read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0)
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

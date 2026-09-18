import os

import matplotlib.pyplot as plt
import numpy as np
import oxbow as ox
import pandas as pd
import pyBigWig
import pypdf
from coolbox.api import *
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer.compute_features_levels import compute_features_levels

from .common import get_cpcdh_intron, merge_adjacent_intervals_with_identical_values
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
            .sort_values(by=["start"], ignore_index=True)
        )

        df_pv = (
            pd
            .concat(
                [
                    pd.DataFrame({
                        "start": [
                            cfg[assemble]["start"],
                            df_pv["end"].to_list()[-1],
                        ],
                        "end": [
                            df_pv["start"].to_list()[0],
                            cfg[assemble]["end"],
                        ],
                        "value": 0.0,
                    }),
                    df_pv,
                    pd.DataFrame({
                        "start": df_pv["end"].to_list()[:-1],
                        "end": df_pv["start"].to_list()[1:],
                        "value": 0.0,
                    }),
                ],
                axis=0,
            )
            .assign(chrom=cfg[assemble]["chrom"])
            .sort_values(by=["start"], ignore_index=True)
        )

        with pyBigWig.open(os.fspath(bw_file), "w") as bw:
            bw.addHeader([(cfg[assemble]["chrom"], cfg[assemble]["length"])])
            bw.addEntries(
                df_pv["chrom"].to_list(),
                df_pv["start"].to_list(),
                ends=df_pv["end"].to_list(),
                values=df_pv["value"].to_list(),
            )


def construct_diff_bw(cfg: dict) -> None:
    df_merge = (
        get_merge_bam(cfg).query("treat.str.endswith('treat')").reset_index(drop=True)
    )

    for exp, protein, treat in zip(
        df_merge["exp"], df_merge["protein"], df_merge["treat"]
    ):
        assemble = treat2assemble(treat)
        chrom = cfg[assemble]["chrom"]
        start = cfg[assemble]["start"]
        end = cfg[assemble]["end"]
        chrom_size = cfg[assemble]["length"]

        treat_bw = cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{treat}.bw"
        control_bw = (
            cfg["data_dir"]
            / "result"
            / "bw"
            / f"{'_'.join(map_to_wild_type_merge(exp, protein, treat))}.bw"
        )
        if not control_bw.exists() or not treat_bw.exists():
            continue

        diff = treat2diff(treat)
        diff_up_bw = cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{diff}.up.bw"
        diff_down_bw = (
            cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{diff}.down.bw"
        )

    with (
        pyBigWig.open(os.fspath(control_bw)) as cb,
        pyBigWig.open(os.fspath(treat_bw)) as tb,
        pyBigWig.open(os.fspath(diff_up_bw), "w") as dub,
        pyBigWig.open(os.fspath(diff_down_bw), "w") as ddb,
    ):
        control_values = cb.values(chrom, start, end, numpy=True)
        treat_values = tb.values(chrom, start, end, numpy=True)
        diff_values = np.nan_to_num(treat_values) - np.nan_to_num(control_values)

        starts, ends, diff_values = merge_adjacent_intervals_with_identical_values(
            starts=np.arange(start, end),
            ends=np.arange(start + 1, end + 1),
            values=diff_values,
        )

        dub.addHeader([(chrom, chrom_size)])
        dub.addEntries(
            [chrom] * len(starts),
            starts,
            ends=ends,
            values=np.maximum(diff_values, 0.0),
        )

        ddb.addHeader([(chrom, chrom_size)])
        ddb.addEntries(
            [chrom] * len(starts),
            starts,
            ends=ends,
            values=-np.minimum(diff_values, 0.0),
        )


def draw_links(
    cfg: dict,
    exp: str,
    protein: str,
    cluster: str,
    assemble: str,
) -> os.PathLike:
    control = "control" if assemble == "hg19" else "mmcontrol"
    if exp != "clip":
        control_f = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_{control}.f.bedpe"
        )
        control_r = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_{control}.r.bedpe"
        )
    else:
        control_f = (
            cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp}_WT_{control}.f.bedpe"
        )
        control_r = (
            cfg["data_dir"] / "result" / "hic" / "bedpe" / f"{exp}_WT_{control}.r.bedpe"
        )

    treat = "delta" if exp != "clip" else "tag"
    if assemble == "mm10":
        treat = f"mm{treat}"

    treat_f = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_{treat}.f.bedpe"
    )
    treat_r = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_{treat}.r.bedpe"
    )

    if (
        not control_f.exists()
        or not control_r.exists()
        or not treat_f.exists()
        or not treat_r.exists()
    ):
        return None

    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)

    chrom = cfg[assemble][cluster]["chrom"]
    start = cfg[assemble][cluster]["start"]
    end = cfg[assemble][cluster]["end"]
    score_to_width = "0.5 + score"  # score is RPM
    diameter_to_height = f"0.5 * max_height * diameter / ({end} - {start})"
    height = 1.5
    tapered = 0.2
    frame = (
        Frame(width=18)
        + XAxis(name=assemble)
        + BEDPE(
            os.fspath(control_f),
            score_to_width=score_to_width,
            diameter_to_height=diameter_to_height,
            tapered=tapered,
            color=cfg["color"]["WT"],
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
            score_to_width=score_to_width,
            diameter_to_height=diameter_to_height,
            tapered=tapered,
            color=cfg["color"]["WT"],
            height=height,
            title="control",
            orientation="inverted",
        )
        + BEDPE(
            os.fspath(treat_f),
            score_to_width=score_to_width,
            diameter_to_height=diameter_to_height,
            tapered=tapered,
            color=cfg["color"][protein],
            height=height,
            title=treat,
        )
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BEDPE(
            os.fspath(treat_r),
            score_to_width=score_to_width,
            diameter_to_height=diameter_to_height,
            tapered=tapered,
            color=cfg["color"][protein],
            height=height,
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
        / f"{assemble}_{exp}_{protein}_{cluster}_links.pdf"
    )

    fig = frame.plot(chrom, start, end)
    fig.savefig(os.fspath(link_file))
    plt.close(fig)

    return link_file


def draw_pre_exons(
    cfg: dict,
    exp: str,
    protein: str,
    cluster: str,
    assemble: str,
) -> os.PathLike:
    control = "control" if assemble == "hg19" else "mmcontrol"
    if exp != "clip":
        control_bw = cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{control}.bw"
    else:
        control_bw = cfg["data_dir"] / "result" / "bw" / f"{exp}_WT_{control}.bw"

    treat = "delta" if exp != "clip" else "tag"
    if assemble == "mm10":
        treat = f"mm{treat}"

    treat_bw = cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{treat}.bw"

    if not control_bw.exists() or not treat_bw.exists():
        return None

    (cfg["data_dir"] / "result" / "hic" / "draw").mkdir(parents=True, exist_ok=True)

    chrom = cfg[assemble][cluster]["chrom"]
    start = cfg[assemble][cluster]["start"]
    end = cfg[assemble][cluster]["end"]

    max_heights = []
    for bwfile in [control_bw, treat_bw]:
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
        + BED(
            os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
            display="collapsed",
            labels=False,
            title=cluster,
        )
        + BigWig(
            os.fspath(control_bw),
            min_value=0,
            max_value=yup,
            color=cfg["color"]["WT"],
            height=height,
            title="control",
        )
        + Spacer(0.5)
        + BigWig(
            os.fspath(treat_bw),
            min_value=0,
            max_value=yup,
            color=cfg["color"][protein],
            height=height,
            title=treat,
        )
        + Spacer(0.5)
        + FrameTitle(protein)
    )
    pre_exon_file = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "draw"
        / f"{assemble}_{exp}_{protein}_{cluster}_pre_exon.pdf"
    )
    fig = frame.plot(chrom, start, end)
    fig.savefig(os.fspath(pre_exon_file))
    plt.close(fig)

    return pre_exon_file


def draw_all(cfg: dict, assemble: str) -> None:
    for exp in ["total", "rna", "clip", "pro"]:
        pdf_files = []
        with pypdf.PdfWriter() as pdf_writer:
            for protein in ["NP220", "MPP8", "PPHLN1", "TASOR"]:
                for cluster in ["alpha"]:
                    pdf_file = draw_links(
                        cfg,
                        exp=exp,
                        protein=protein,
                        cluster=cluster,
                        assemble=assemble,
                    )
                    if pdf_file is not None:
                        pdf_writer.append(pdf_file)
                        pdf_files.append(pdf_file)

                    pdf_file = draw_pre_exons(
                        cfg,
                        exp=exp,
                        protein=protein,
                        cluster=cluster,
                        assemble=assemble,
                    )
                    if pdf_file is not None:
                        pdf_writer.append(pdf_file)
                        pdf_files.append(pdf_file)

            if pdf_files:
                pdf_writer.write(
                    cfg["data_dir"]
                    / "result"
                    / "hic"
                    / "draw"
                    / f"{assemble}_{exp}.pdf"
                )

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
    assemble: str,
):
    control = "control" if assemble == "hg19" else "mmcontrol"
    if exp != "clip":
        control_f = (
            cfg["data_dir"]
            / "bam"
            / "merge"
            / "precursor"
            / f"{exp}_{protein}_{control}.f.bam"
        )
        control_r = (
            cfg["data_dir"]
            / "bam"
            / "merge"
            / "precursor"
            / f"{exp}_{protein}_{control}.r.bam"
        )
    else:
        control_f = (
            cfg["data_dir"]
            / "bam"
            / "merge"
            / "precursor"
            / f"{exp}_WT_{control}.f.bam"
        )
        control_r = (
            cfg["data_dir"]
            / "bam"
            / "merge"
            / "precursor"
            / f"{exp}_WT_{control}.r.bam"
        )

    treat = "delta" if exp != "clip" else "tag"
    if assemble == "mm10":
        treat = f"mm{treat}"

    treat_f = (
        cfg["data_dir"]
        / "bam"
        / "merge"
        / "precursor"
        / f"{exp}_{protein}_{treat}.f.bam"
    )
    treat_r = (
        cfg["data_dir"]
        / "bam"
        / "merge"
        / "precursor"
        / f"{exp}_{protein}_{treat}.r.bam"
    )

    if (
        not control_f.exists()
        or not control_r.exists()
        or not treat_f.exists()
        or not treat_r.exists()
    ):
        return

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

    for name, se, pos in zip(df_se["name"], df_se["se"], df_se["pos"]):
        chrom = cfg[assemble]["chrom"]
        start = pos - cfg["read_length"]
        end = pos + cfg["read_length"]
        frame = (
            XAxis(name=assemble)
            + BAM(
                os.fspath(control_f),
                length_ratio_thresh=cfg["length_ratio_thresh"],
                color=cfg["color"]["WT"],
                height=estimate_height(control_f, chrom, start, end),
                title="control",
            )
            + BED(
                os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
                display="collapsed",
                labels=False,
                title=f"{name}:{se}",
            )
            + BAM(
                os.fspath(control_r),
                length_ratio_thresh=cfg["length_ratio_thresh"],
                color=cfg["color"]["WT"],
                height=estimate_height(control_r, chrom, start, end),
                title="control",
                orientation="inverted",
            )
            + BAM(
                os.fspath(treat_f),
                length_ratio_thresh=cfg["length_ratio_thresh"],
                color=cfg["color"][protein],
                height=estimate_height(treat_f, chrom, start, end),
                title=treat,
            )
            + BED(
                os.fspath(cfg["data_dir"] / "result" / f"{assemble}.12.bed"),
                display="collapsed",
                labels=False,
                title=f"{name}:{se}",
            )
            + BAM(
                os.fspath(treat_r),
                length_ratio_thresh=cfg["length_ratio_thresh"],
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


def draw_reads_all(cfg: dict, assemble: str):
    for exp in ["total", "rna", "pro", "clip"]:
        pdf_files = []
        with pypdf.PdfWriter() as pdf_writer:
            for protein in ["NP220", "MPP8", "PPHLN1", "TASOR"]:
                for pdf_file in draw_reads(cfg, exp, protein, assemble):
                    pdf_writer.append(pdf_file)
                    pdf_files.append(pdf_file)

            if pdf_files:
                pdf_writer.write(
                    cfg["data_dir"]
                    / "result"
                    / "hic"
                    / "draw"
                    / f"{assemble}_{exp}_reads.pdf"
                )

        for pdf_file in pdf_files:
            pdf_file.unlink()

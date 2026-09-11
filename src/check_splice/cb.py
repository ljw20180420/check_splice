import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import oxbow as ox
import pandas as pd
import py2bit
import pyBigWig
import pypdf
import pysam
from coolbox.api import *
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer.compute_features_levels import compute_features_levels

from .utils import (
    bw_merge_adjacent_intervals_with_identical_values,
    get_precursor_pos,
    select_total_count,
)


def pairs_to_bedpe(cfg: dict) -> None:
    shutil.rmtree(cfg["data_dir"] / "result" / "hic" / "bedpe", ignore_errors=True)
    (cfg["data_dir"] / "result" / "hic" / "bedpe").mkdir(parents=True, exist_ok=True)
    for pairs_file in os.listdir(cfg["data_dir"] / "result" / "hic" / "pairs"):
        exp, protein, treat = pairs_file.removesuffix(".pairs").split("_", 3)
        total_count = select_total_count(cfg, exp, protein, treat)

        pairs_file = cfg["data_dir"] / "result" / "hic" / "pairs" / pairs_file

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

        assemble = "mm10" if treat.startswith("mm") else "hg19"
        with py2bit.open(cfg[assemble]["2bit"]) as tb:
            donors = []
            acceptors = []
            for chrom1, start1, chrom2, start2 in zip(
                df["chrom1"], df["start1"], df["chrom2"], df["start2"]
            ):
                assert chrom1 == chrom2 and start1 < start2, "illegal order"
                donors.append(tb.sequence(chrom1, start1, start1 + 2))
                acceptors.append(tb.sequence(chrom2, start2 - 2, start2))

        df = df.assign(
            donor=donors,
            acceptor=acceptors,
        )

        bedpe_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / pairs_file.with_suffix(".bedpe").name
        )

        df.query("donor == 'GT' and acceptor == 'AG'").drop(
            columns=["donor", "acceptor"]
        ).to_csv(
            bedpe_file.with_suffix(".f.bedpe"), sep="\t", header=False, index=False
        )

        df.query("donor == 'CT' and acceptor == 'AC'").drop(
            columns=["donor", "acceptor"]
        ).to_csv(
            bedpe_file.with_suffix(".r.bedpe"), sep="\t", header=False, index=False
        )

        df.query(
            "(donor != 'GT' or acceptor != 'AG') and (donor != 'CT' or acceptor != 'AC')"
        ).drop(columns=["donor", "acceptor"]).to_csv(
            bedpe_file.with_suffix(".o.bedpe"), sep="\t", header=False, index=False
        )


def diff_bedpe(cfg: dict, exp: str, protein: str, orientation: str) -> None:
    treat = "delta" if exp != "clip" else "tag"
    treat_file = (
        cfg["data_dir"]
        / "result"
        / "hic"
        / "bedpe"
        / f"{exp}_{protein}_{treat}.{orientation}.bedpe"
    )
    if exp != "clip":
        control_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_{protein}_control.{orientation}.bedpe"
        )
    else:
        control_file = (
            cfg["data_dir"]
            / "result"
            / "hic"
            / "bedpe"
            / f"{exp}_WT_control.{orientation}.bedpe"
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
            "donor",
            "acceptor",
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
            "donor",
            "acceptor",
        ],
    )
    df = df_treat.merge(
        df_control,
        on=[
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "strand1",
            "strand2",
        ],
        how="outer",
    ).assign(
        name=lambda df: df["name_x"].combine_first(df["name_y"]),
        score=lambda df: df["score_x"].fillna(0) - df["score_y"].fillna(0),
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
            / f"{exp}_{protein}_diff.up.{orientation}.bedpe"
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
            / f"{exp}_{protein}_diff.down.{orientation}.bedpe"
        ),
        sep="\t",
        index=False,
        header=False,
    )


def diff_bedpe_all(cfg: dict) -> None:
    for exp in ["total", "rna", "pro", "clip"]:
        for protein in ["NP220", "MPP8", "PPHLN1", "TASOR"]:
            for orientation in ["f", "r"]:
                diff_bedpe(cfg, exp, protein, orientation)


def get_exon_pre(cfg: dict, assemble: str):
    assert assemble in ["hg19", "mm10"], "unknown assemble"
    df_cpcdh = (
        pd
        .read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0)
        .query("type=='exon' and name.str.lower().str.startswith('pcdha')")
        .reset_index(drop=True)
    )
    df_cpcdh = df_cpcdh.assign(
        pre_start=lambda df: (
            [df.loc[0, "start"].item() - cfg["size_before_first"]]
            + df["end"].to_list()[:-1]
        ),
        pre_end=lambda df: df["start"],
    )

    dfs = []
    for exp in ["total", "rna", "pro", "clip"]:
        for protein in ["WT", "NP220", "MPP8", "PPHLN1", "TASOR"]:
            for wt in [True, False]:
                if wt:
                    treat = "control"
                else:
                    if exp != "clip":
                        treat = "delta"
                    else:
                        treat = "tag"

                if assemble == "mm10":
                    treat = f"mm{treat}"

                bam_file = (
                    cfg["data_dir"] / "bam" / "merge" / f"{exp}_{protein}_{treat}.bam"
                )

                if not bam_file.exists():
                    continue

                with pysam.AlignmentFile(os.fspath(bam_file)) as sam:
                    exon_counts = []
                    pre_counts = []
                    for chrom, start, end, pre_start, pre_end in zip(
                        df_cpcdh["chrom"],
                        df_cpcdh["start"],
                        df_cpcdh["end"],
                        df_cpcdh["pre_start"],
                        df_cpcdh["pre_end"],
                    ):
                        exon_count = 0
                        for read in sam.fetch(chrom, start, end):
                            if read.is_secondary:
                                continue
                            if not read.is_mapped:
                                continue
                            if read.is_supplementary:
                                continue

                            exon_count += 1

                        exon_counts.append(exon_count)

                        pre_count = 0
                        for read in sam.fetch(chrom, pre_start, pre_end):
                            if read.is_secondary:
                                continue
                            if not read.is_mapped:
                                continue
                            if read.is_supplementary:
                                continue

                            pre_count += 1

                        pre_counts.append(pre_count)

                dfs.append(
                    df_cpcdh.copy()[
                        ["chrom", "start", "end", "name", "pre_start", "pre_end"]
                    ].assign(
                        exon_count=exon_counts,
                        pre_count=pre_counts,
                        exp=exp,
                        protein=protein,
                        treat=treat,
                    )
                )

    pd.concat(dfs, ignore_index=True).to_csv(
        cfg["data_dir"] / "result" / f"{assemble}_exon_pre.csv", index=False
    )


def construct_artifact_bw(cfg: dict, assemble: str) -> None:
    df_splice = (
        pd
        .read_csv(cfg["data_dir"] / "result" / f"{assemble}_splice.csv", header=0)
        .query("name.str.lower().str.startswith('pcdha')")
        .reset_index(drop=True)
    )
    for exp in ["total", "rna", "pro", "clip"]:
        for protein in ["WT", "NP220", "MPP8", "PPHLN1", "TASOR"]:
            for wt in [True, False]:
                if wt:
                    treat = "control"
                else:
                    if exp != "clip":
                        treat = "delta"
                    else:
                        treat = "tag"

                if assemble == "mm10":
                    treat = f"mm{treat}"

                df_splice_slice = df_splice.query(
                    "exp == @exp and protein == @protein and treat == @treat"
                ).assign(**{
                    "splice %": lambda df: (
                        df["connect"] / (df["connect"] + df["cover.start"]) * 100
                    )
                })

                if len(df_splice_slice) == 0:
                    continue

                (cfg["data_dir"] / "result" / "bw").mkdir(exist_ok=True, parents=True)
                bw_file = (
                    cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{treat}.bw"
                )

                df_pv = (
                    pd
                    .DataFrame({
                        "start": df_splice_slice["start"].to_list(),
                        "value": df_splice_slice["splice %"].fillna(0.0).to_list(),
                    })
                    .assign(
                        end=lambda df: df["start"] + 1,
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
                    .sort_values(by=["chrom", "start"], ignore_index=True)
                )

                with pyBigWig.open(os.fspath(bw_file), "w") as bw:
                    bw.addHeader([(cfg[assemble]["chrom"], cfg[assemble]["length"])])
                    bw.addEntries(
                        df_pv["chrom"].to_list(),
                        df_pv["start"].to_list(),
                        ends=df_pv["end"].to_list(),
                        values=df_pv["value"].to_list(),
                    )


def construct_diff_bw(cfg: dict, assemble: str) -> None:
    assert assemble in ["hg19", "mm10"], "unknown assemble"
    control = "control" if assemble == "hg19" else "mmcontrol"
    chrom = cfg[assemble]["chrom"]
    start = cfg[assemble]["start"]
    end = cfg[assemble]["end"]
    for exp in ["total", "rna", "pro", "clip"]:
        for protein in ["WT", "NP220", "MPP8", "PPHLN1", "TASOR"]:
            if exp != "clip":
                control_bw = (
                    cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{control}.bw"
                )
            else:
                control_bw = (
                    cfg["data_dir"] / "result" / "bw" / f"{exp}_WT_{control}.bw"
                )

            treat = "delta" if exp != "clip" else "tag"
            if assemble == "mm10":
                treat = f"mm{treat}"

            treat_bw = cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{treat}.bw"

            if not control_bw.exists() or not treat_bw.exists():
                continue

            diff = "diff" if assemble == "hg19" else "mmdiff"
            diff_up_bw = (
                cfg["data_dir"] / "result" / "bw" / f"{exp}_{protein}_{diff}.up.bw"
            )
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
                diff_values = np.nan_to_num(treat_values) - np.nan_to_num(
                    control_values
                )

                starts, ends, diff_values = (
                    bw_merge_adjacent_intervals_with_identical_values(
                        starts=np.arange(start, end),
                        ends=np.arange(start + 1, end + 1),
                        values=diff_values,
                    )
                )

                dub.addHeader([(cfg[assemble]["chrom"], cfg[assemble]["length"])])
                dub.addEntries(
                    [cfg[assemble]["chrom"]] * len(starts),
                    starts,
                    ends=ends,
                    values=np.maximum(diff_values, 0.0),
                )

                ddb.addHeader([(cfg[assemble]["chrom"], cfg[assemble]["length"])])
                ddb.addEntries(
                    [cfg[assemble]["chrom"]] * len(starts),
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
        get_precursor_pos(cfg, assemble)
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

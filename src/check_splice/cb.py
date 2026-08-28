import os

import matplotlib.pyplot as plt
import pandas as pd
import pypdf
from coolbox.api import *


def get_bedpe_file_name(
    cfg: dict, exp_protein_wt: str, orientations: str
) -> os.PathLike:
    if orientations == ["ff"] or orientations == ["rr"]:
        orientation = orientations[0]
    elif orientations == ["ff", "rr"] or orientations == ["rr", "ff"]:
        orientation = "fr"
    else:
        raise ValueError("Invalid orientations")

    bedpe_file = f"{orientation}.bedpe"

    if exp_protein_wt:
        bedpe_file = f"{exp_protein_wt}_{bedpe_file}"

    return cfg["data_dir"] / "result" / "hic" / bedpe_file, orientation


def draw_links(
    cfg: dict,
    exp_protein_wt: str,
    orientations: str,
):
    dfs = []
    for orientation in orientations:
        if exp_protein_wt:
            pairs_file = (
                cfg["data_dir"]
                / "result"
                / "hic"
                / f"{exp_protein_wt}_{orientation}.pairs"
            )
        else:
            pairs_file = cfg["data_dir"] / "result" / "hic" / f"{orientation}.pairs"

        dfs.append(
            pd.read_csv(
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
        )

    bedpe_file, orientation = get_bedpe_file_name(cfg, exp_protein_wt, orientations)

    df = pd.concat(dfs, ignore_index=True)
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
        .agg(name=pd.NamedAgg("readID", "first"), score=pd.NamedAgg("readID", "count"))
        .reset_index()[
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

    regions = {
        "cpcdh": {
            "chrom": cfg["chrom"],
            "start": cfg["start"],
            "end": cfg["end"],
        },
        "alpha": cfg["alpha"],
        "beta": cfg["beta"],
        "gamma": cfg["gamma"],
    }
    for cluster, region in regions.items():
        frame = (
            XAxis(name="hg19")
            + BEDPE(os.fspath(bedpe_file))
            + TrackHeight(5)
            + Title("splice")
            + Spacer(1)
            + BED(
                os.fspath(cfg["data_dir"] / "result" / "hg19.12.bed"),
                display="interlaced",
                labels=False,
            )
            + Title("gene")
            + Spacer(1)
            + BED(
                os.fspath(cfg["data_dir"] / "result" / "pCBS.bed"),
                display="interlaced",
                labels=False,
            )
            + Title("pCBS")
            + FrameTitle(f"{exp_protein_wt}:{orientation}:{cluster}")
        )
        fig = frame.plot(f"{region['chrom']}:{region['start']}-{region['end']}")
        fig.savefig(os.fspath(bedpe_file.with_suffix(f".{cluster}.pdf")))
        plt.close(fig)

        yield bedpe_file.with_suffix(f".{cluster}.pdf")


def draw_links_all(cfg):
    pdf_files = []
    with pypdf.PdfWriter() as pdf_writer:
        for orientations in [["ff", "rr"], ["ff"], ["rr"]]:
            for pdf_file in draw_links(
                cfg,
                "",
                orientations,
            ):
                pdf_writer.append(pdf_file)
                pdf_files.append(pdf_file)
            with open("exp_protein_wts.txt", "r") as fd:
                for exp_protein_wt in fd:
                    exp_protein_wt = exp_protein_wt.strip()
                    for pdf_file in draw_links(cfg, exp_protein_wt, orientations):
                        pdf_writer.append(pdf_file)
                        pdf_files.append(pdf_file)

        pdf_writer.write(cfg["data_dir"] / "result" / "hic" / "links.pdf")

    for pdf_file in pdf_files:
        pdf_file.unlink()

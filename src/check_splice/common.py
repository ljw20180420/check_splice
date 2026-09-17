import os
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pysam
import sh


def get_bam_read_count(bamfile: os.PathLike) -> int:
    with pysam.AlignmentFile(os.fspath(bamfile)) as bam:
        total_count = sum(
            1
            for read in bam
            if not read.is_secondary and read.is_mapped and not read.is_supplementary
        )

    return total_count


def prepare_gene_bed12(
    gtffile: os.PathLike, outfile: os.PathLike, addtional_filter: Callable
) -> None:
    gtffile = Path(os.fspath(gtffile))
    outfile = Path(os.fspath(outfile))

    gtfToGenePred = sh.Command("gtfToGenePred")
    genePredToBigGenePred = sh.Command("genePredToBigGenePred")
    gtfToGenePred(
        "-genePredExt",
        os.fspath(gtffile),
        os.fspath(gtffile.with_suffix(".gp")),
    )
    df_gp = pd.read_csv(gtffile.with_suffix(".gp"), sep="\t", header=None)
    df_gp[[11] + list(range(1, 15))].to_csv(
        gtffile.with_suffix(".gp"),
        sep="\t",
        header=False,
        index=False,
    )
    genePredToBigGenePred(
        os.fspath(gtffile.with_suffix(".gp")),
        os.fspath(gtffile.with_suffix(".bgp")),
    )
    df_bgp = (
        pd
        .read_csv(gtffile.with_suffix(".bgp"), sep="\t", header=None)[list(range(12))]
        .rename(
            columns={
                0: "chrom",
                1: "chromStart",
                2: "chromEnd",
                3: "name",
                4: "score",
                5: "strand",
                6: "thickStart",
                7: "thickEnd",
                8: "itemRgb",
                9: "blockCount",
                10: "blockSizes",
                11: "blockStarts",
            }
        )
        .sort_values(by=["chrom", "chromStart"], ignore_index=True)
    )

    df_bgp = (
        df_bgp
        .query("not name.str.lower().str.startswith('pcdha') or blockCount == 4")
        .query(
            "not name.str.lower().str.startswith('pcdhb') or blockCount == 1 or name.str.lower() == 'pcdhb9'"
        )
        .query("not name.str.lower().str.startswith('pcdhb@')")
        .query("not name.str.lower().str.startswith('pcdhg') or blockCount == 4")
        .reset_index(drop=True)
    )

    df_bgp = addtional_filter(df_bgp)

    df_bgp.to_csv(
        outfile,
        sep="\t",
        header=False,
        index=False,
    )
    outfile.with_suffix(".bed.bgz").unlink(missing_ok=True)
    outfile.with_suffix(".bed.bgz.tbi").unlink(missing_ok=True)


def get_cpcdh_exon(gtffile: os.PathLike, chrom: str) -> pd.DataFrame:
    df = pd.read_csv(
        gtffile,
        sep="\t",
        names=[
            "chrom",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ],
    )
    df = df.query(
        "chrom == @chrom and attributes.str.lower().str.contains(r'pcdh[abg][abc]?[0-9]{1,2}') and (feature=='exon' or feature=='CDS')"
    ).reset_index(drop=True)
    attributes = df["attributes"].str.split(expand=True)
    df = df.assign(
        start=lambda df: df["start"] - 1,
        name=attributes[9].str.strip('";'),
        transcript_id=attributes[3].str.strip('";'),
        exon_number=attributes[5].str.strip('";').astype(int),
        total_exon_number=lambda df: df.groupby("transcript_id")[
            "exon_number"
        ].transform(max),
    ).drop(
        columns=[
            "source",
            "score",
            "frame",
            "attributes",
        ]
    )
    df = (
        df
        .query(
            "(exon_number == 1 and (total_exon_number == 4 or name.str.lower().str.startswith('pcdhb'))) or (exon_number > 1 and (name.str.lower() == 'pcdha1' or name.str.lower() == 'pcdhga1'))"
        )
        .reset_index(drop=True)
        .assign(
            repeat=lambda df: df.groupby(["feature", "name", "exon_number"])[
                "name"
            ].transform("count"),
        )
        .query("repeat == 1 or transcript_id.str.contains(r'(?:^NM_018|NM_002)')")
        .reset_index(drop=True)
        .assign(
            name=lambda df: df.apply(
                lambda row: (
                    row["name"]
                    if row["exon_number"] == 1
                    else f"ace{row['exon_number'] - 1}"
                    if row["name"].lower() == "pcdha1"
                    else f"gce{row['exon_number'] - 1}"
                ),
                axis=1,
            )
        )
    )
    df = df.pivot_table(
        values=["start", "end"],
        index=["chrom", "strand", "name", "transcript_id"],
        columns="feature",
    )
    df.columns = df.columns.to_flat_index().map(lambda tp: f"{tp[1]}_{tp[0]}")
    df = (
        df
        .reset_index()
        .rename(columns={"exon_start": "start", "exon_end": "end"})
        .assign(
            CDS_start=lambda df: df["CDS_start"].fillna("."),
            CDS_end=lambda df: df["CDS_end"].fillna("."),
        )
        .assign(score=".")[
            [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "CDS_start",
                "CDS_end",
                "transcript_id",
            ]
        ]
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["start", "end"], ignore_index=True)
    )

    return df


def get_cpcdh_intron(cpcdh_csv: os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(cpcdh_csv, header=0)

    intron_names = []
    intron_starts = []
    intron_ends = []

    intron_end = df.query("name == 'ace1'")["start"].item()
    for name, intron_start in df.query("name.str.lower().str.startswith('pcdha')")[
        ["name", "end"]
    ].itertuples(index=False):
        intron_names.append(f"{name}_ace1")
        intron_starts.append(intron_start)
        intron_ends.append(intron_end)

    intron_names.append("ace1_ace2")
    intron_starts.append(df.query("name == 'ace1'")["end"].item())
    intron_ends.append(df.query("name == 'ace2'")["start"].item())

    intron_names.append("ace2_ace3")
    intron_starts.append(df.query("name == 'ace2'")["end"].item())
    intron_ends.append(df.query("name == 'ace3'")["start"].item())

    intron_end = df.query("name == 'gce1'")["start"].item()
    for name, intron_start in df.query("name.str.lower().str.startswith('pcdhg')")[
        ["name", "end"]
    ].itertuples(index=False):
        intron_names.append(f"{name}_gce1")
        intron_starts.append(intron_start)
        intron_ends.append(intron_end)

    intron_names.append("gce1_gce2")
    intron_starts.append(df.query("name == 'gce1'")["end"].item())
    intron_ends.append(df.query("name == 'gce2'")["start"].item())

    intron_names.append("gce2_gce3")
    intron_starts.append(df.query("name == 'gce2'")["end"].item())
    intron_ends.append(df.query("name == 'gce3'")["start"].item())

    return (
        pd
        .DataFrame({
            "start": intron_starts,
            "end": intron_ends,
            "name": intron_names,
        })
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["start"], ignore_index=True)
    )

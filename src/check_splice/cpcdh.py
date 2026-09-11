import os

import pandas as pd
import sh


def prepare_gene_bed12(cfg: dict) -> None:
    gtfToGenePred = sh.Command("gtfToGenePred")
    genePredToBigGenePred = sh.Command("genePredToBigGenePred")
    for gtf_file in ["hg19.ncbiRefSeq.gtf", "mm10.ncbiRefSeq.gtf"]:
        gtf_file = cfg["data_dir"] / "data" / gtf_file
        gtfToGenePred(
            "-genePredExt",
            os.fspath(gtf_file),
            os.fspath(gtf_file.with_suffix(".gp")),
        )
        df_gp = pd.read_csv(gtf_file.with_suffix(".gp"), sep="\t", header=None)
        df_gp[[11] + list(range(1, 15))].to_csv(
            gtf_file.with_suffix(".gp"),
            sep="\t",
            header=False,
            index=False,
        )
        genePredToBigGenePred(
            os.fspath(gtf_file.with_suffix(".gp")),
            os.fspath(gtf_file.with_suffix(".bgp")),
        )
        df_bgp = (
            pd
            .read_csv(gtf_file.with_suffix(".bgp"), sep="\t", header=None)[
                list(range(12))
            ]
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

        if gtf_file.name == "hg19.ncbiRefSeq.gtf":
            df_bgp = (
                df_bgp
                .query("not name.str.startswith('PCDHA') or blockCount == 4")
                .query(
                    "not name.str.startswith('PCDHB') or blockCount == 1 or name == 'PCDHB9'"
                )
                .query("not name.str.startswith('PCDHB@')")
                .query("not name.str.startswith('PCDHG') or blockCount == 4")
                .query("name != 'PCDHA1' or blockSizes.str.startswith('2545')")
                .query("name != 'PCDHA6' or blockSizes.str.startswith('2526')")
                .query("name != 'PCDHA10' or blockSizes.str.startswith('2540')")
                .query("name != 'PCDHGA11' or blockSizes.str.startswith('2610')")
                .query("name != 'PCDHGC3' or blockSizes.str.startswith('2581')")
                .query(
                    "name != 'LOC112267934' and name != 'LOC101926905' and name != 'LOC100419552' and name != 'SLC25A2' and name != 'TAF7' and name != 'RN7SL68P'"
                )
                .reset_index(drop=True)
            )
        elif gtf_file.name == "mm10.ncbiRefSeq.gtf":
            df_bgp = (
                df_bgp
                .query("not name.str.startswith('Pcdha') or blockCount == 4")
                .query("not name.str.startswith('Pcdhb') or blockCount == 1")
                .query("not name.str.startswith('PCDHG') or blockCount == 4")
                .query(
                    "name != 'Gm37013' and name != 'Gm38666' and name != 'Gm38667' and name != 'Gm36858' and name != 'Gm18150' and name != 'Gm19035' and name != 'Gm18529' and name != 'Gm20162' and name != 'Slc25a2' and name != 'Taf7' and name != 'Gm8242' and name != 'Gm24401' and name != 'Gm29994'"
                )
                .reset_index(drop=True)
            )
        else:
            raise ValueError("unknown gtf")

        out_file = (
            "hg19.12.bed" if gtf_file.name == "hg19.ncbiRefSeq.gtf" else "mm10.12.bed"
        )
        out_file = cfg["data_dir"] / "result" / out_file

        df_bgp.to_csv(
            out_file,
            sep="\t",
            header=False,
            index=False,
        )
        out_file.with_suffix(".bed.bgz").unlink(missing_ok=True)
        out_file.with_suffix(".bed.bgz.tbi").unlink(missing_ok=True)


def get_hg19_cpcdh_exon(gtffile: os.PathLike) -> pd.DataFrame:
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
        "chrom == 'chr5' and attributes.str.contains(r'PCDH[ABG][ABC]?[0-9]{1,2}') and (feature=='exon' or feature=='CDS')"
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
            "exon_number == 1 and (total_exon_number == 4 or name.str.contains(r'^PCDHB'))"
        )
        .reset_index(drop=True)
        .assign(
            repeat=lambda df: df.groupby(["feature", "name"])["name"].transform(
                "count"
            ),
        )
        .query("repeat == 1 or transcript_id.str.contains(r'(?:^NM_018|NM_002)')")
        .reset_index(drop=True)
        .drop(columns=["exon_number", "total_exon_number"])
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
    )
    starts = [
        140358533,
        140362059,
        140389211,
        140874373,
        140884959,
        140890513,
    ]
    ends = [
        140358592,
        140362148,
        140391932,
        140874432,
        140885048,
        140892542,
    ]
    df = (
        pd
        .concat(
            [
                df,
                pd.DataFrame({
                    "chrom": ["chr5"] * 6,
                    "start": starts,
                    "end": ends,
                    "name": ["ace1", "ace2", "ace3", "gce1", "gce2", "gce3"],
                    "score": ["."] * 6,
                    "strand": ["+"] * 6,
                    "CDS_start": starts,
                    "CDS_end": ends,
                    "transcript_id": ["ace1", "ace2", "ace3", "gce1", "gce2", "gce3"],
                }),
            ],
        )
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["start", "end"], ignore_index=True)
    )

    return df


def get_mm10_cpcdh_exon(gtffile: os.PathLike) -> pd.DataFrame:
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
        "chrom == 'chr18' and attributes.str.contains(r'Pcdh[abg][abc]?[0-9]{1,2}') and (feature=='exon' or feature=='CDS')"
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
            "(exon_number == 1 and (total_exon_number == 4 or name.str.startswith('Pcdhb'))) or (exon_number > 1 and (name == 'Pcdha1' or name == 'Pcdhga1'))"
        )
        .reset_index(drop=True)
        .assign(
            repeat=lambda df: df.groupby(["feature", "name", "exon_number"])[
                "name"
            ].transform("count"),
        )
        .query("repeat == 1")
        .reset_index(drop=True)
        .assign(
            name=lambda df: df.apply(
                lambda row: (
                    row["name"]
                    if row["exon_number"] == 1
                    else f"ace{row['exon_number'] - 1}"
                    if row["name"] == "Pcdha1"
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


def get_cpcdh_intron(df: pd.DataFrame) -> pd.DataFrame:
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

    df["type"] = "exon"
    chrom = df.loc[0, "chrom"]
    df_intron = pd.DataFrame({
        "chrom": chrom,
        "start": intron_starts,
        "end": intron_ends,
        "name": intron_names,
        "score": ".",
        "strand": "+",
        "CDS_start": ".",
        "CDS_end": ".",
        "transcript_id": ".",
        "type": "intron",
    })

    return (
        pd
        .concat([df, df_intron], ignore_index=True)
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["type", "start", "end"], ignore_index=True)
    )


def get_cpcdh(cfg: dict) -> None:
    df = get_hg19_cpcdh_exon(cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gtf.gz")
    df = get_cpcdh_intron(df)
    (cfg["data_dir"] / "result").mkdir(exist_ok=True, parents=True)
    df.to_csv(cfg["data_dir"] / "result" / "hg19_cpcdh.csv", index=False)

    df = get_mm10_cpcdh_exon(cfg["data_dir"] / "data" / "mm10.ncbiRefSeq.gtf.gz")
    df = get_cpcdh_intron(df)
    (cfg["data_dir"] / "result").mkdir(exist_ok=True, parents=True)
    df.to_csv(cfg["data_dir"] / "result" / "mm10_cpcdh.csv", index=False)

import os

import pandas as pd


def get_cpcdh_exon(gtffile: os.PathLike) -> pd.DataFrame:
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
        .astype({
            "start": "int64",
            "end": "int64",
        })
        .sort_values(by=["start", "end"], ignore_index=True)
    )

    return df


def get_cpcdh_intron(df: pd.DataFrame) -> pd.DataFrame:
    intron_names = []
    intron_starts = []
    intron_ends = []

    intron_end = df.query("name == 'ace1'")["start"].item()
    for name, intron_start in df.query("name.str.startswith('PCDHA')")[
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
    for name, intron_start in df.query("name.str.startswith('PCDHG')")[
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
    df_intron = pd.DataFrame({
        "chrom": "chr5",
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

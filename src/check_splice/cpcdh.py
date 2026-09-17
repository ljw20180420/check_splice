import pandas as pd


def gene_bed12_hg19_cpcdh_filter(df_bgp: pd.DataFrame) -> pd.DataFrame:
    return (
        df_bgp
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


def gene_bed12_mm10_cpcdh_filter(df_bgp: pd.DataFrame) -> pd.DataFrame:
    return df_bgp.query(
        "name != 'Gm37013' and name != 'Gm38666' and name != 'Gm38667' and name != 'Gm36858' and name != 'Gm18150' and name != 'Gm19035' and name != 'Gm18529' and name != 'Gm20162' and name != 'Slc25a2' and name != 'Taf7' and name != 'Gm8242' and name != 'Gm24401' and name != 'Gm29994'"
    ).reset_index(drop=True)


def fix_mm10_Pcdhgb8_CDS(df: pd.DataFrame) -> pd.DataFrame:
    df.loc[df["name"] == "Pcdhgb8", "CDS_start"] = 37761878
    df.loc[df["name"] == "Pcdhgb8", "CDS_end"] = 37764299

    return df

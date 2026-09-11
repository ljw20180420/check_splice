import pandas as pd


def get_pCBS(cfg: dict) -> None:
    for assemble in ["hg19", "mm10"]:
        shift_file = f"{assemble}_pCBS_shift.csv"
        cpcdh_file = cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"

        df_shift = pd.read_csv(shift_file, header=0)
        df_cpcdh = pd.read_csv(cpcdh_file, header=0)
        df = df_shift.merge(
            right=df_cpcdh,
            how="inner",
            left_on="gene",
            right_on="name",
            validate="1:1",
        )

        df = df.astype({"CDS_start": float, "CDS_end": float}).astype({
            "CDS_start": int,
            "CDS_end": int,
        })
        df = df.assign(
            end=lambda df: df["CDS_start"] + df["shift"] + 2,
            start=lambda df: df["end"] - 27,
        )[["chrom", "start", "end", "name"]].assign(score=".", strand="+")

        out_file = cfg["data_dir"] / "result" / f"{assemble}_pCBS.bed"
        df.to_csv(
            out_file,
            sep="\t",
            header=False,
            index=False,
        )
        out_file.with_suffix(".bed.bgz").unlink(missing_ok=True)
        out_file.with_suffix(".bed.bgz.tbi").unlink(missing_ok=True)

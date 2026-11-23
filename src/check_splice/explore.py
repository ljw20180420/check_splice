import numpy as np
import pandas as pd
import plotly.express as px
import py2bit
from Bio.Seq import Seq


def check_RS_at_exon_tail(cfg: dict):
    for assemble in ["hg19", "mm10"]:
        df = pd.read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0)
        with py2bit.open(cfg[assemble]["2bit"]) as tb:
            for chrom, end, name in zip(df["chrom"], df["end"], df["name"]):
                seq1 = tb.sequence(chrom, end - 10, end)
                seq2 = tb.sequence(chrom, end, end + 10)
                print(seq1, seq2, assemble, name)


def recursive_microhomology(seq1: str, seq2: str) -> pd.DataFrame:
    mask_matrix = np.pad(
        array=np.array(list(seq1))[:, None] == np.array(list(seq2)),
        pad_width=((1, 0), (1, 0)),
        mode="constant",
        constant_values=False,
    )
    mh_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1))
    for i in range(1, len(mh_matrix)):
        current_mask = mask_matrix[i]
        shifted_mask = np.append(current_mask[1:], False)
        mask_matrix[i, current_mask] = mask_matrix[i - 1, shifted_mask] + 1
        mask_matrix[i - 1, shifted_mask] = 0

    rows, cols = np.nonzero(mh_matrix)

    return pd.DataFrame({
        "row": rows,
        "col": cols,
        "mh": mh_matrix[rows, cols],
    })


def check_microhomology_of_two_strands_of_cpcdh_exons(cfg: dict):
    for assemble in ["hg19", "mm10"]:
        df = pd.read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0)
        with py2bit.open(cfg[assemble]["2bit"]) as tb:
            for chrom, start, end, name in zip(
                df["chrom"], df["start"], df["end"], df["name"]
            ):
                seq1 = tb.sequence(chrom, start, end)
                seq2 = Seq(seq1).reverse_complement()

                df = recursive_microhomology(seq1, seq2)
                xs = []
                ys = []
                for row, col, mh in df.itertuples(index=False):
                    xs.extend([col - mh, col, None])
                    ys.extend([row - mh, row, None])

                fig = px.line(
                    pd.DataFrame({
                        "x": xs,
                        "y": ys,
                    }),
                    x="x",
                    y="y",
                    labels={
                        "x": "crick",
                        "y": "watson",
                    },
                )

                (cfg["data_dir"] / "explore").mkdir(exist_ok=True, parents=True)
                fig.write_image(cfg["data_dir"] / "explore" / f"{name}.pdf")

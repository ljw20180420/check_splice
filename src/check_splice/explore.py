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
    seq1 = seq1.upper()
    seq2 = seq2.upper()
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
        mh_matrix[i, current_mask] = mh_matrix[i - 1, shifted_mask] + 1
        mh_matrix[i - 1, shifted_mask] = 0

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
                seq1 = tb.sequence(chrom, start, end).upper()
                seq2 = str(Seq(seq1).reverse_complement())

                df_mh = recursive_microhomology(seq1, seq2)
                df_mh = df_mh.query("mh >= @cfg['explore']['min_mh']").reset_index(
                    drop=True
                )
                xs = []
                ys = []
                for row, col, mh in df_mh.itertuples(index=False):
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
                    title=f"{assemble}:{name}",
                    width=800,
                    height=800,
                )
                fig.update_xaxes(
                    side="top",
                    range=[0, len(seq2)],
                    constrain="domain",
                )
                fig.update_yaxes(
                    scaleanchor="x",
                    scaleratio=1,
                    range=[len(seq1), 0],
                    constrain="domain",
                )
                fig.update_traces(line_width=0.5)

                (cfg["data_dir"] / "explore").mkdir(exist_ok=True, parents=True)
                fig.write_image(cfg["data_dir"] / "explore" / f"{assemble}_{name}.pdf")

                yield cfg["data_dir"] / "explore" / f"{assemble}_{name}.pdf"

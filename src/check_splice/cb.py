import os

import matplotlib.pyplot as plt
from coolbox.api import *


def draw_links(cfg: dict, pairs_file: os.PathLike):
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
            XAxis()
            + Pairs(os.fspath(pairs_file))
            + Title("splice")
            + BED(os.fspath(cfg["data_dir"] / "data" / "hg19.12.bed"))
            + Title("gene")
            + BED(os.fspath(cfg["data_dir"] / "result" / "pCBS.bed"))
            + Title("pCBS")
        )
        fig = frame.plot(f"{region['chrom']}:{region['start']}-{region['end']}")
        fig.savefig(os.fspath(pairs_file.with_suffix(f".{cluster}.pdf")))
        plt.close(fig)

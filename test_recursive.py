#!/usr/bin/env python

import pandas as pd
import py2bit

from check_splice import config

cfg = config.pcdh()

for assemble in ["hg19", "mm10"]:
    df = pd.read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0)
    with py2bit.open(cfg[assemble]["2bit"]) as tb:
        for chrom, end, name in zip(df["chrom"], df["end"], df["name"]):
            seq1 = tb.sequence(chrom, end - 10, end)
            seq2 = tb.sequence(chrom, end, end + 10)
            print(seq1, seq2, assemble, name)

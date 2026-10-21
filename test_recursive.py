#!/usr/bin/env python

import pandas as pd
import py2bit

from check_splice import config

cfg = config.pcdh()

for assemble in ["hg19", "mm10"]:
    df = pd.read_csv(cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv", header=0)
    with py2bit.open(cfg[assemble]["2bit"]) as tb:
        for chrom, start, end, name in zip(
            df["chrom"], df["start"], df["end"], df["name"]
        ):
            seq1 = tb.sequence(chrom, start - 10, start + 2)
            seq2 = tb.sequence(chrom, end - 2, end + 5)
            print(seq1, seq2, assemble, name)

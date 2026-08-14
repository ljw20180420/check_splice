#!/usr/bin/env python

from check_splice import config, cpcdh

cfg = config.pcdh()
df = cpcdh.get_cpcdh_exon(cfg["data_dir"] / "hg19.ncbiRefSeq.gtf.gz")
df = cpcdh.get_cpcdh_intron(df)
df.to_csv(cfg["data_dir"] / "cpcdh.csv", index=False)

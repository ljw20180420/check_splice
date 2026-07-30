#!/usr/bin/env python

from check_splice import config, get_cpcdh

cfg = config.pcdh()
df = get_cpcdh(cfg["data_dir"] / "hg19.ncbiRefSeq.gtf.gz")
df.to_csv(cfg["data_dir"] / "cpcdh.csv", index=False)

#!/usr/bin/env python

from check_splice import config, pcbs

cfg = config.pcdh()
df = pcbs.get_pCBS(
    shift_file=cfg["data_dir"] / "data" / "pCBS_shift.csv",
    cpcdh_file=cfg["data_dir"] / "result" / "cpcdh.csv",
)
df.to_csv(cfg["data_dir"] / "result" / "pCBS.bed", sep="\t", header=False, index=False)

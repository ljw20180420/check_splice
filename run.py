#!/usr/bin/env python

from check_splice import get_cpcdh

df = get_cpcdh("hg19.ncbiRefSeq.gtf.gz")
df.to_csv("/home/ljw/sdc1/hush/cpcdh.csv", index=False)

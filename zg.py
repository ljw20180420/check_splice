#!/usr/bin/env python

import sys

import py2bit
from Bio.Seq import Seq

from check_splice import config

chrom = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
strand = sys.argv[4]
assemble = sys.argv[5]


cfg = config.pcdh()
with py2bit.open(cfg[assemble]["2bit"]) as tb:
    seq = tb.sequence(chrom, start, end)
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())

print(seq)

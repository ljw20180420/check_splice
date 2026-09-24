#!/usr/bin/env python

import io
import sys

import pandas as pd
import py2bit
import sh
from Bio.Seq import Seq

chrom1 = sys.argv[1]
start1 = int(sys.argv[2])
end1 = int(sys.argv[3])
strand1 = sys.argv[4]
chrom2 = sys.argv[5]
start2 = int(sys.argv[6])
end2 = int(sys.argv[7])
strand2 = sys.argv[8]

db = "/home/ljw/sdb1/ucsc/hubs/myHub/lmm10/lmm10.2bit"
with py2bit.open(db) as tb:
    seq1 = tb.sequence(chrom1, start1, end1)
    if strand1 == "-":
        seq1 = str(Seq(seq1).reverse_complement())
    seq2 = tb.sequence(chrom2, start2, end2)
    if strand2 == "-":
        seq2 = str(Seq(seq2).reverse_complement())

    seq = seq1 + seq2

print(len(seq))
blat = sh.Command("blat")
result = blat(
    "-out=blast8",
    "-stepSize=5",
    "-repMatch=2253",
    "-minScore=0",
    "-minIdentity=0",
    f"{db}:chr18:36923470-37904446",
    "stdin",
    "stdout",
    _in=f">seq\n{seq}\n",
)

df = pd.read_csv(
    io.StringIO(result),
    sep="\t",
    names=[
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ],
)

print(df)

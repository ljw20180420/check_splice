#!/usr/bin/env python

import pandas as pd

from check_splice import config, count_splice

cfg = config.pcdh()
cpcdh_pair = pd.read_csv(cfg["data_dir"] / "cpcdh_pair.csv")
for nearby_file in [
    "total_rna_seq/WT.bam.chimeric.nearby",
    "rna_seq/RNA-seq-PPdel-PPHLN1.merge.bam.chimeric.nearby",
    "rna_seq/RNA-seq-WT-PPHLN1.merge.bam.chimeric.nearby",
    "rna_seq/TASORdel.bam.chimeric.nearby",
]:
    df_nearby = pd.read_csv(cfg["data_dir"] / nearby_file, header=0)
    counts = []
    for seqname1, end1, strand1, seqname2, start2, strand2 in zip(
        cpcdh_pair["seqname1"],
        cpcdh_pair["end1"],
        cpcdh_pair["strand1"],
        cpcdh_pair["seqname2"],
        cpcdh_pair["start2"],
        cpcdh_pair["strand2"],
    ):
        count = count_splice(
            df_nearby,
            splice1_chrom=seqname1,
            splice1_pos=end1,
            splice1_strand=strand1,
            splice2_chrom=seqname2,
            splice2_pos=start2,
            splice2_strand=strand2,
            splice_thres=cfg["splice_thres"],
        )
        counts.append(count)

    cpcdh_pair.assign(
        count=counts,
        size1=lambda df: df["end1"] - df["start1"],
        size2=lambda df: df["end2"] - df["start2"],
    ).to_csv(cfg["data_dir"] / f"{nearby_file}.splice", index=False)

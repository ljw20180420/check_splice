#!/usr/bin/env python

import pandas as pd

from check_splice import config, count_non_splice

cfg = config.pcdh()
cpcdh = pd.read_csv(cfg["data_dir"] / "cpcdh.csv")
for nearby_file in [
    "total_rna_seq/WT.bam.chimeric.nearby",
    "rna_seq/RNA-seq-PPdel-PPHLN1.merge.bam.chimeric.nearby",
    "rna_seq/RNA-seq-WT-PPHLN1.merge.bam.chimeric.nearby",
    "rna_seq/TASORdel.bam.chimeric.nearby",
]:
    df_nearby = pd.read_csv(cfg["data_dir"] / nearby_file, header=0)
    counts = []
    for seqname, end, strand in zip(
        cpcdh["seqname"],
        cpcdh["end"],
        cpcdh["strand"],
    ):
        count = count_non_splice(
            df_nearby,
            splice_chrom=seqname,
            splice_pos=end,
            splice_strand=strand,
            non_splice_thres=cfg["non_splice_thres"],
        )
        counts.append(count)

    cpcdh.assign(count=counts).query("count > 0").to_csv(
        cfg["data_dir"] / f"{nearby_file}.non_splice", index=False
    )

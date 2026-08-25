from pathlib import Path


def pcdh():
    return {
        "chrom": "chr5",
        "start": 140125000,
        "end": 140922000,
        "alpha": {
            "chrom": "chr5",
            "start": 140125000,
            "end": 140414000,
        },
        "beta": {
            "chrom": "chr5",
            "start": 140414000,
            "end": 140668000,
        },
        "gamma": {
            "chrom": "chr5",
            "start": 140668000,
            "end": 140922000,
        },
        "data_dir": Path("/home/ljw/sdc1/hush"),
        "cover_threshold": 3,
        "tss_extend": 30,
        "exon_end_extend": 30,
    }

from pathlib import Path


def pcdh():
    return {
        "chrom": "chr5",
        "start": 140125000,
        "end": 140922000,
        "data_dir": Path("/home/ljw/sdc1/hush"),
        "splice_thres": 0,
        "non_splice_thres": 3,
    }

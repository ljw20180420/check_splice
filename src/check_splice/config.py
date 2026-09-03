from pathlib import Path


def pcdh():
    return {
        "color": {
            "WT": "#000000",
            "NP220": "#FF1493",
            "PPHLN1": "#A52A2A",
            "MPP8": "#800080",
            "TASOR": "#008000",
            "INCREASE": "#FF0000",
            "DECREASE": "#0000FF",
        },
        "flip": "R1",  # which reads to flip (R1 or R2)
        "chrom": "chr5",
        "start": 140158536,
        "end": 140964431,
        "alpha": {
            "chrom": "chr5",
            "start": 140158536,
            "end": 140425184,
        },
        "beta": {
            "chrom": "chr5",
            "start": 140425185,
            "end": 140703220,
        },
        "gamma": {
            "chrom": "chr5",
            "start": 140703221,
            "end": 140964431,
        },
        "data_dir": Path("/home/ljw/sdc1/hush"),
        "cover_threshold": 3,
        "tss_extend": 30,
        "exon_end_extend": 30,
    }

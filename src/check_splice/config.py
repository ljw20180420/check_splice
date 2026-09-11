from pathlib import Path


def pcdh():
    return {
        "size_before_first": 3000,
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
        "hg19": {
            "chrom": "chr5",
            "length": 180915260,
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
            "2bit": "/home/ljw/sdb1/ucsc/hubs/myHub/lhg19/lhg19.2bit",
        },
        "mm10": {
            "chrom": "chr18",
            "length": 90702639,
            "start": 36927858,
            "end": 37842465,
            "alpha": {
                "chrom": "chr18",
                "start": 36927858,
                "end": 37194516,
            },
            "beta": {
                "chrom": "chr18",
                "start": 37194516,
                "end": 37540035,
            },
            "gamma": {
                "chrom": "chr18",
                "start": 37540035,
                "end": 37842465,
            },
            "2bit": "/home/ljw/sdb1/ucsc/hubs/myHub/lmm10/lmm10.2bit",
        },
        "data_dir": Path("/home/ljw/sdc1/hush"),
        "cover_threshold": 6,
        "tss_extend": 30,
        "exon_end_extend": 30,
    }

from pathlib import Path


def pcdh():
    return {
        "length_ratio_thresh": 1e-5,
        "read_length": 150,
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
            "start": 36923470,
            "end": 37904446,
            "alpha": {
                "chrom": "chr18",
                "start": 36923470,
                "end": 37232199,
            },
            "beta": {
                "chrom": "chr18",
                "start": 37232200,
                "end": 37653385,
            },
            "gamma": {
                "chrom": "chr18",
                "start": 37653386,
                "end": 37904446,
            },
            "2bit": "/home/ljw/sdb1/ucsc/hubs/myHub/lmm10/lmm10.2bit",
        },
        "data_dir": Path("/home/ljw/sdc1/hush"),
        "cover_threshold": 6,
        "tss_extend": 30,
        "exon_end_extend": 30,
    }

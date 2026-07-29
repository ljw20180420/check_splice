import os

import pandas as pd
import pysam


def is_overlap(
    chrom1: str, start1: int, end1: int, chrom2: str, start2: int, end2: int
) -> bool:
    if chrom1 != chrom2:
        return False
    return min(end1, end2) > max(start1, start2)


def parse_sa(read: pysam.AlignedSegment):
    sa_tag = read.get_tag("SA")
    for alignment_str in sa_tag.split(";"):
        if not alignment_str:
            continue

        chrom, start, strand, cigar, _ = alignment_str.split(",")
        start = int(start) - 1
        temp_read = pysam.AlignedSegment()
        temp_read.cigarstring = cigar
        end = start + temp_read.reference_length

        yield chrom, start, end, strand


def parse_chimeric(
    samfile: os.PathLike, locus_chrom: str, locus_start: int, locus_end: int
) -> pd.DataFrame:
    with pysam.AlignmentFile(os.fspath(samfile)) as fd:
        json_lines = []
        for read in fd.fetch():
            if read.is_secondary:
                continue
            if not read.is_mapped:
                continue
            if read.is_supplementary:
                continue

            segs = []

            seg_chrom, seg_start, seg_end = (
                read.reference_name,
                read.reference_start,
                read.reference_end,
            )
            seg_strand = "+" if read.is_forward else "-"
            if not is_overlap(
                seg_chrom,
                seg_start,
                seg_end,
                locus_chrom,
                locus_start,
                locus_end,
            ):
                continue
            segs.append(f"{seg_chrom}:{seg_start}:{seg_end}:{seg_strand}")

            if read.has_tag("SA"):
                overlap = True
                for seg_chrom, seg_start, seg_end, seg_strand in parse_sa(read):
                    if not is_overlap(
                        seg_chrom,
                        seg_start,
                        seg_end,
                        locus_chrom,
                        locus_start,
                        locus_end,
                    ):
                        overlap = False
                        break
                    segs.append(f"{seg_chrom}:{seg_start}:{seg_end}:{seg_strand}")

                if not overlap:
                    continue

            json_lines.append({
                "query_name": read.query_name,
                "is_read1": read.is_read1,
                "is_qcfail": read.is_qcfail,
                "is_duplicate": read.is_duplicate,
                "mapping_quality": read.mapping_quality,
                "segs": ";".join(segs),
            })

    return pd.DataFrame(json_lines)


def get_cpcdh(gtffile: os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(
        gtffile,
        sep="\t",
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ],
    )
    df = df.query(
        "seqname == 'chr5' and attributes.str.contains(r'PCDH[ABG][ABC]?[0-9]{1,2}') and feature=='exon' and attributes.str.contains('exon_number \"1\"') and end - start >= 1500"
    ).reset_index(drop=True)
    attributes = df["attributes"].str.split(expand=True)
    df = df.assign(
        gene_name=attributes[9].str.strip('";'),
        transcript_id=attributes[3].str.strip('";'),
    ).drop(columns=["attributes"])
    return df


# def nearby_explode(df: pd.DataFrame) -> pd.DataFrame:

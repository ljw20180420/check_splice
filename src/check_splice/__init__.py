import os

import pandas as pd
import pysam


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
        "seqname == 'chr5' and attributes.str.contains(r'PCDH[ABG][ABC]?[0-9]{1,2}') and feature=='exon'"
    ).reset_index(drop=True)
    attributes = df["attributes"].str.split(expand=True)
    df = df.assign(
        gene_name=attributes[9].str.strip('";'),
        transcript_id=attributes[3].str.strip('";'),
        exon_number=attributes[5].str.strip('";').astype(int),
        total_exon_number=lambda df: df.groupby("transcript_id")[
            "exon_number"
        ].transform(max),
    ).drop(
        columns=[
            "source",
            "feature",
            "score",
            "frame",
            "attributes",
        ]
    )
    df = (
        df
        .query(
            "exon_number == 1 and (total_exon_number == 4 or gene_name.str.contains(r'^PCDHB'))"
        )
        .reset_index(drop=True)
        .assign(
            repeat=lambda df: df.groupby("gene_name")["gene_name"].transform("count"),
        )
        .query("repeat == 1 or transcript_id.str.contains(r'(?:^NM_018|NM_002)')")
        .reset_index(drop=True)
        .drop(columns=["exon_number", "total_exon_number", "repeat"])
    )
    df = pd.concat(
        [
            df,
            pd.DataFrame({
                "seqname": ["chr5"] * 6,
                "start": [
                    140358533,
                    140362059,
                    140389211,
                    140874373,
                    140884959,
                    140890513,
                ],
                "end": [
                    140358592,
                    140362148,
                    140391932,
                    140874432,
                    140885048,
                    140892542,
                ],
                "strand": ["+"] * 6,
                "gene_name": ["ace1", "ace2", "ace3", "gce1", "gce2", "gce3"],
                "transcript_id": ["ace1", "ace2", "ace3", "gce1", "gce2", "gce3"],
            }),
        ],
    ).sort_values(by=["start", "end"], ignore_index=True)

    return df


def pair_cpcdh(cpcdh_file: os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(cpcdh_file)
    df = (
        df
        .merge(df, how="cross", suffixes=["1", "2"])
        .query("gene_name1 != gene_name2 and end1 < start2")
        .reset_index(drop=True)
    )

    return df


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


def get_intron_only_blocks(read):
    """
    Returns a list of (start, end) tuples split strictly by 'N' (introns).
    Deletions (D) and insertions (I) are treated as continuous segments.
    """
    # 0-based genomic starting position of the alignment
    current_pos = read.reference_start
    blocks = []

    # Track the start of the current block
    block_start = current_pos

    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

    for op, length in read.cigartuples:
        if op in (0, 2, 7, 8):  # Operators that consume reference genome space
            current_pos += length

        elif op == 3:  # Intron / Reference Skip (N)
            # End the current block before the intron starts
            if current_pos > block_start:
                blocks.append((block_start, current_pos))

            # Skip past the intron region
            current_pos += length
            # Set the start of the next block to the end of the intron
            block_start = current_pos

        elif op not in (
            1,
            4,
            5,
        ):  # Insertions and clips do not move the reference cursor
            raise ValueError("unknown cigar operation")

    # Append the final block after parsing the last CIGAR operation
    if current_pos > block_start:
        blocks.append((block_start, current_pos))

    return blocks


def parse_chimeric(
    samfile: os.PathLike, locus_chrom: str, locus_start: int, locus_end: int
) -> pd.DataFrame:
    with pysam.AlignmentFile(os.fspath(samfile)) as fd:
        json_lines = []
        for read in fd.fetch(
            contig=locus_chrom,
            start=locus_start,
            end=locus_end,
        ):
            if read.is_secondary:
                continue
            if not read.is_mapped:
                continue
            if read.is_supplementary:
                continue

            segs = []
            seg_chrom = read.reference_name
            seg_strand = "+" if read.is_forward else "-"

            blocks = get_intron_only_blocks(read)
            for seg_start, seg_end in blocks:
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
                "is_chimeric": len(segs) > 1,
                "has_SA": read.has_tag("SA"),
            })

    return pd.DataFrame(json_lines)


def nearby_explode(chimeric_file: os.PathLike) -> pd.DataFrame:
    def nearby(segs):
        if len(segs) == 1:
            return [f"{segs[0]}:{segs[0]}"]
        else:
            return [f"{segs[i]}:{segs[i + 1]}" for i in range(len(segs) - 1)]

    df = pd.read_csv(chimeric_file, header=0)
    df = (
        df
        .assign(
            segs=lambda df: df["segs"].str.split(";"),
            nearby_segs=lambda df: df["segs"].map(nearby),
        )
        .drop(columns=["segs"])
        .explode(column="nearby_segs", ignore_index=True)
    )
    nearby_segs = df["nearby_segs"].str.split(":", expand=True)
    df = df.assign(
        chrom1=nearby_segs[0],
        start1=nearby_segs[1],
        end1=nearby_segs[2],
        strand1=nearby_segs[3],
        chrom2=nearby_segs[4],
        start2=nearby_segs[5],
        end2=nearby_segs[6],
        strand2=nearby_segs[7],
    ).drop(columns=["nearby_segs"])

    return df


def count_splice(
    df_nearby: pd.DataFrame,
    splice1_chrom: str,
    splice1_pos: int,
    splice1_strand: str,
    splice2_chrom: str,
    splice2_pos: int,
    splice2_strand,
    splice_thres: int,
) -> int:
    df = df_nearby.query("is_chimeric").reset_index(drop=True)
    df = df.assign(
        splice1=lambda df: (
            (df["chrom1"] == splice1_chrom)
            & (df["strand1"] == splice1_strand)
            & ((df["end1"] - splice1_pos).abs() <= splice_thres)
        ),
        splice2=lambda df: (
            (df["chrom2"] == splice2_chrom)
            & (df["strand2"] == splice2_strand)
            & ((df["start2"] - splice2_pos).abs() <= splice_thres)
        ),
        splice=lambda df: df["splice1"] & df["splice2"],
    )
    return df["splice"].sum().item()


def count_non_splice(
    df_cpcdh: pd.DataFrame,
    splice_chrom: str,
    splice_pos: int,
    splice_strand: str,
    non_splice_thres: int,
) -> int:
    df = df_cpcdh.assign(
        cover1=lambda df: (
            (df["chrom1"] == splice_chrom)
            & (df["strand1"] == splice_strand)
            & (splice_pos >= df["start1"] + non_splice_thres)
            & (splice_pos <= df["end1"] - non_splice_thres)
        ),
        cover2=lambda df: (
            (df["chrom2"] == splice_chrom)
            & (df["strand2"] == splice_strand)
            & (splice_pos >= df["start2"] + non_splice_thres)
            & (splice_pos <= df["end2"] - non_splice_thres)
        ),
        cover=lambda df: df["cover1"] | df["cover2"],
    )
    df = df.query("cover").drop_duplicates(subset=["query_name", "is_read1"])

    return len(df)

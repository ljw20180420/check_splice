import io
import os
import re
from collections.abc import Callable, Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import pypdf
import pysam
import sh
from matplotlib import cm, colors


def get_bam_read_count(bamfile: os.PathLike) -> int:
    with pysam.AlignmentFile(os.fspath(bamfile)) as bam:
        total_count = sum(
            1
            for read in bam
            if not read.is_secondary and read.is_mapped and not read.is_supplementary
        )

    return total_count


def prepare_gene_bed12(
    gtffile: os.PathLike, outfile: os.PathLike, addtional_filter: Callable
) -> None:
    gtffile = Path(os.fspath(gtffile))
    outfile = Path(os.fspath(outfile))

    gtfToGenePred = sh.Command("gtfToGenePred")
    genePredToBigGenePred = sh.Command("genePredToBigGenePred")
    gtfToGenePred(
        "-genePredExt",
        os.fspath(gtffile),
        os.fspath(gtffile.with_suffix(".gp")),
    )
    df_gp = pd.read_csv(gtffile.with_suffix(".gp"), sep="\t", header=None)
    df_gp[[11] + list(range(1, 15))].to_csv(
        gtffile.with_suffix(".gp"),
        sep="\t",
        header=False,
        index=False,
    )
    genePredToBigGenePred(
        os.fspath(gtffile.with_suffix(".gp")),
        os.fspath(gtffile.with_suffix(".bgp")),
    )
    df_bgp = (
        pd
        .read_csv(gtffile.with_suffix(".bgp"), sep="\t", header=None)[list(range(12))]
        .rename(
            columns={
                0: "chrom",
                1: "chromStart",
                2: "chromEnd",
                3: "name",
                4: "score",
                5: "strand",
                6: "thickStart",
                7: "thickEnd",
                8: "itemRgb",
                9: "blockCount",
                10: "blockSizes",
                11: "blockStarts",
            }
        )
        .sort_values(by=["chrom", "chromStart"], ignore_index=True)
    )

    df_bgp = (
        df_bgp
        .query("not name.str.lower().str.startswith('pcdha') or blockCount == 4")
        .query(
            "not name.str.lower().str.startswith('pcdhb') or blockCount == 1 or name.str.lower() == 'pcdhb9'"
        )
        .query("not name.str.lower().str.startswith('pcdhb@')")
        .query("not name.str.lower().str.startswith('pcdhg') or blockCount == 4")
        .reset_index(drop=True)
    )

    df_bgp = addtional_filter(df_bgp)

    df_bgp.to_csv(
        outfile,
        sep="\t",
        header=False,
        index=False,
    )
    outfile.with_suffix(".bed.bgz").unlink(missing_ok=True)
    outfile.with_suffix(".bed.bgz.tbi").unlink(missing_ok=True)


def get_cpcdh_exon(gtffile: os.PathLike, chrom: str) -> pd.DataFrame:
    df = pd.read_csv(
        gtffile,
        sep="\t",
        names=[
            "chrom",
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
        "chrom == @chrom and attributes.str.lower().str.contains(r'pcdh[abg][abc]?[0-9]{1,2}') and (feature=='exon' or feature=='CDS')"
    ).reset_index(drop=True)
    attributes = df["attributes"].str.split(expand=True)
    df = df.assign(
        start=lambda df: df["start"] - 1,
        name=attributes[9].str.strip('";'),
        transcript_id=attributes[3].str.strip('";'),
        exon_number=attributes[5].str.strip('";').astype(int),
        total_exon_number=lambda df: df.groupby("transcript_id")[
            "exon_number"
        ].transform(max),
    ).drop(
        columns=[
            "source",
            "score",
            "frame",
            "attributes",
        ]
    )
    df = (
        df
        .query(
            "(exon_number == 1 and (total_exon_number == 4 or name.str.lower().str.startswith('pcdhb'))) or (exon_number > 1 and (name.str.lower() == 'pcdha1' or name.str.lower() == 'pcdhga1'))"
        )
        .reset_index(drop=True)
        .assign(
            repeat=lambda df: df.groupby(["feature", "name", "exon_number"])[
                "name"
            ].transform("count"),
        )
        .query("repeat == 1 or transcript_id.str.contains(r'(?:^NM_018|NM_002)')")
        .reset_index(drop=True)
        .assign(
            name=lambda df: df.apply(
                lambda row: (
                    row["name"]
                    if row["exon_number"] == 1
                    else f"ace{row['exon_number'] - 1}"
                    if row["name"].lower() == "pcdha1"
                    else f"gce{row['exon_number'] - 1}"
                ),
                axis=1,
            )
        )
    )
    df = df.pivot_table(
        values=["start", "end"],
        index=["chrom", "strand", "name", "transcript_id"],
        columns="feature",
    )
    df.columns = df.columns.to_flat_index().map(lambda tp: f"{tp[1]}_{tp[0]}")
    df = (
        df
        .reset_index()
        .rename(columns={"exon_start": "start", "exon_end": "end"})
        .assign(
            CDS_start=lambda df: df["CDS_start"].fillna("."),
            CDS_end=lambda df: df["CDS_end"].fillna("."),
        )
        .assign(score=".")[
            [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "CDS_start",
                "CDS_end",
                "transcript_id",
            ]
        ]
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["start", "end"], ignore_index=True)
    )

    return df


def get_cpcdh_intron(cpcdh_csv: os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(cpcdh_csv, header=0)

    intron_names = []
    intron_starts = []
    intron_ends = []

    intron_end = df.query("name == 'ace1'")["start"].item()
    for name, intron_start in df.query("name.str.lower().str.startswith('pcdha')")[
        ["name", "end"]
    ].itertuples(index=False):
        intron_names.append(f"{name}_ace1")
        intron_starts.append(intron_start)
        intron_ends.append(intron_end)

    intron_names.append("ace1_ace2")
    intron_starts.append(df.query("name == 'ace1'")["end"].item())
    intron_ends.append(df.query("name == 'ace2'")["start"].item())

    intron_names.append("ace2_ace3")
    intron_starts.append(df.query("name == 'ace2'")["end"].item())
    intron_ends.append(df.query("name == 'ace3'")["start"].item())

    intron_end = df.query("name == 'gce1'")["start"].item()
    for name, intron_start in df.query("name.str.lower().str.startswith('pcdhg')")[
        ["name", "end"]
    ].itertuples(index=False):
        intron_names.append(f"{name}_gce1")
        intron_starts.append(intron_start)
        intron_ends.append(intron_end)

    intron_names.append("gce1_gce2")
    intron_starts.append(df.query("name == 'gce1'")["end"].item())
    intron_ends.append(df.query("name == 'gce2'")["start"].item())

    intron_names.append("gce2_gce3")
    intron_starts.append(df.query("name == 'gce2'")["end"].item())
    intron_ends.append(df.query("name == 'gce3'")["start"].item())

    return (
        pd
        .DataFrame({
            "chrom": df.loc[0, "chrom"],
            "start": intron_starts,
            "end": intron_ends,
            "name": intron_names,
        })
        .astype({"start": "int64", "end": "int64"})
        .sort_values(by=["start"], ignore_index=True)
    )


def is_mapped_primary_first(read: pysam.AlignedSegment) -> bool:
    if read.is_secondary:
        return False
    if not read.is_mapped:
        return False

    return not read.is_supplementary


def filter_bam_reads(bamfile: os.PathLike, chrom: str, start: int, end: int):
    with pysam.AlignmentFile(os.fspath(bamfile)) as fd:
        for read in fd.fetch(
            contig=chrom,
            start=start,
            end=end,
        ):
            if is_mapped_primary_first(read):
                yield read


class ParseSamRead:
    def __init__(self):
        self.cigar_parser = re.compile(r"(\d+)([MIDNSHP=XB])")
        self.align_sting_parser = re.compile(r"\d+[MID]")

    def flip_to_query_raw_strand(
        self,
        query_block_start: int,
        query_block_end: int,
        align_string: list[str],
        strand: str,
        query_length: int,
    ):
        if strand == "+":
            return (
                query_block_start,
                query_block_end,
                "".join(align_string),
            )
        else:
            # strand == "-"
            return (
                query_length - query_block_end,
                query_length - query_block_start,
                "".join(reversed(align_string)),
            )

    def parse_cigar(self, start: int, cigarstring: str, strand: str, query_length: int):
        ref_current_pos = start
        ref_block_start = start
        query_current_pos = 0
        length_ops = [
            (int(length), op) for length, op in self.cigar_parser.findall(cigarstring)
        ]
        if length_ops[0][1] == "S":
            query_current_pos += length_ops[0][0]
            length_ops = length_ops[1:]
        query_block_start = query_current_pos

        align_string = []
        for length, op in length_ops:
            if op in ("=", "X"):
                op = "M"

            if op in ("M", "I", "D"):
                align_string.append(f"{length}{op}")
                if op in ("M", "D"):
                    ref_current_pos += length
                if op in ("M", "I"):
                    query_current_pos += length

            elif op == "N":
                yield (
                    ref_block_start,
                    ref_current_pos,
                    *self.flip_to_query_raw_strand(
                        query_block_start,
                        query_current_pos,
                        align_string,
                        strand,
                        query_length,
                    ),
                )

                ref_current_pos += length
                ref_block_start = ref_current_pos
                query_block_start = query_current_pos
                align_string = []

        yield (
            ref_block_start,
            ref_current_pos,
            *self.flip_to_query_raw_strand(
                query_block_start,
                query_current_pos,
                align_string,
                strand,
                query_length,
            ),
        )

    def parse_sa(self, sa_tag: str):
        for alignment_str in sa_tag.split(";"):
            if not alignment_str:
                continue

            chrom, start, strand, cigar, _ = alignment_str.split(",")
            start = int(start) - 1

            yield chrom, start, strand, cigar

    def parse_read(self, read: pysam.AlignedSegment):
        chroms = [read.reference_name]
        starts = [read.reference_start]
        strands = ["+" if read.is_forward else "-"]
        cigars = [read.cigarstring]
        if read.has_tag("SA"):
            for chrom, start, strand, cigar in self.parse_sa(read.get_tag("SA")):
                chroms.append(chrom)
                starts.append(start)
                strands.append(strand)
                cigars.append(cigar)

        for chrom, start, strand, cigar in zip(chroms, starts, strands, cigars):
            for (
                ref_block_start,
                ref_block_end,
                query_block_start,
                query_block_end,
                align_string,
            ) in self.parse_cigar(start, cigar, strand, read.query_length):
                yield (
                    chrom,
                    ref_block_start,
                    ref_block_end,
                    strand,
                    query_block_start,
                    query_block_end,
                    align_string,
                )

    def flip_read(
        self,
        ref_block_chrom: str,
        ref_block_start: int,
        ref_block_end: int,
        ref_block_strand: str,
        query_block_start: int,
        query_block_end: int,
        align_string: str,
        query_length: int,
    ):
        return (
            ref_block_chrom,
            ref_block_start,
            ref_block_end,
            "+" if ref_block_strand == "-" else "-",
            query_length - query_block_end,
            query_length - query_block_start,
            "".join(reversed(self.align_sting_parser.findall(align_string))),
        )


def pair_to_hic(
    pair_file: os.PathLike, resolutions: list[int], chrom_sizes: os.PathLike
) -> None:
    hic_file = pair_file.with_suffix(".hic")
    hictk = sh.Command("hictk")
    hictk(
        "load",
        "--format",
        "4dn",
        "--bin-size",
        f"{resolutions[0]}",
        "--chrom-sizes",
        os.fspath(chrom_sizes),
        "--force",
        os.fspath(pair_file),
        os.fspath(hic_file),
    )

    hictk(
        "zoomify",
        "--resolutions",
        *[f"{resolution}" for resolution in resolutions],
        "--force",
        os.fspath(hic_file),
        f"{os.fspath(hic_file.with_suffix('.m.hic'))}",
    )

    sh.mv(os.fspath(hic_file.with_suffix(".m.hic")), os.fspath(hic_file))

    hictk("balance", "scale", os.fspath(hic_file))


def get_bw(
    bamfile: os.PathLike,
    total_count: int,
    bin_size: int,
    chrom: str,
    start: int,
    end: int,
    chrom_size,
) -> None:
    bamfile = Path(bamfile)
    bamCoverage = sh.Command("bamCoverage")
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        mapped_count = bam.mapped
    if mapped_count > 0:
        bamCoverage(
            "--bam",
            os.fspath(bamfile),
            "-o",
            os.fspath(bamfile.with_suffix(".bw")),
            "-r",
            f"{chrom}:{start}:{end}",
            "--binSize",
            bin_size,
            "--scaleFactor",
            1_000_000 / total_count,
        )
    else:
        mid = (start + end) // 2
        with pyBigWig.open(os.fspath(bamfile.with_suffix(".bw")), "w") as bw:
            bw.addHeader([(chrom, chrom_size)])
            bw.addEntries(chrom, [mid], values=[0.0], span=1)


def merge_adjacent_intervals_with_identical_values(
    starts: np.ndarray, ends: np.ndarray, values: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    change_indices = np.where(np.diff(values) != 0)[0] + 1

    ends = np.concatenate((starts[change_indices], [ends[-1]]))
    starts = starts[np.concatenate(([0], change_indices))]
    values = values[np.concatenate(([0], change_indices))]

    return starts, ends, values


def substract_bigwig(
    bigwig_file1: os.PathLike,
    bigwig_file2: os.PathLike,
    chrom: str,
    start: int,
    end: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with (
        pyBigWig.open(os.fspath(bigwig_file1)) as bw1,
        pyBigWig.open(os.fspath(bigwig_file2)) as bw2,
    ):
        bw1_values = bw1.values(chrom, start, end, numpy=True)
        bw2_values = bw2.values(chrom, start, end, numpy=True)
        diff_values = np.nan_to_num(bw1_values) - np.nan_to_num(bw2_values)

        starts, ends, diff_values = merge_adjacent_intervals_with_identical_values(
            starts=np.arange(start, end),
            ends=np.arange(start + 1, end + 1),
            values=diff_values,
        )

    df = pd.DataFrame({"start": starts, "end": ends, "diff_value": diff_values}).query(
        "diff_value != 0"
    )

    return df["start"].to_numpy(), df["end"].to_numpy(), df["diff_value"].to_numpy()


def write_bigwig(
    chrom: str,
    chrom_size: int,
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    bigwig_file: os.PathLike,
):
    if len(starts) == 0:
        starts = np.append(starts, [0, 1])
        ends = np.append(ends, [1, 2])
        values = np.append(values, [0, 0])

    with pyBigWig.open(os.fspath(bigwig_file), "w") as bw:
        bw.addHeader([(chrom, chrom_size)])
        bw.addEntries(
            [chrom] * len(starts),
            starts,
            ends=ends,
            values=values,
        )


def read_bedpe(bedpe_file: os.PathLike) -> pd.DataFrame:
    return pd.read_csv(
        bedpe_file,
        sep="\t",
        names=[
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "name",
            "score",
            "strand1",
            "strand2",
        ],
    )


def substract_bedpe(bedpe_file1: os.PathLike, bedpe_file2: os.PathLike) -> pd.DataFrame:
    df1 = read_bedpe(bedpe_file1)
    df2 = read_bedpe(bedpe_file2)

    df = df1.merge(
        df2,
        on=[
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "strand1",
            "strand2",
        ],
        how="outer",
    )

    df = df.assign(
        name=lambda df: df["name_x"].fillna("") + "-" + df["name_y"].fillna(""),
        score=lambda df: df["score_x"].fillna(0) - df["score_y"].fillna(0),
    )[
        [
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "name",
            "score",
            "strand1",
            "strand2",
        ]
    ]

    return df


def summation_bedpe(
    bedpe_files: list[os.PathLike], total_counts: list[int]
) -> pd.DataFrame:
    dfs = [
        read_bedpe(bedpe_file).assign(
            score=lambda df, total_count=total_count: df["score"] * total_count
        )
        for bedpe_file, total_count in zip(bedpe_files, total_counts)
    ]

    df_sum = dfs[0]
    for df in dfs[1:]:
        df_sum = df_sum.merge(
            df,
            on=[
                "chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "strand1",
                "strand2",
            ],
            how="outer",
        )

        df_sum = df_sum.assign(
            name=lambda df: df["name_x"].fillna("") + ":" + df["name_y"].fillna(""),
            score=lambda df: df["score_x"].fillna(0) + df["score_y"].fillna(0),
        )
        df_sum = df_sum.assign(
            name=lambda df: df["name"].str.strip(":"),
        )[
            [
                "chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "name",
                "score",
                "strand1",
                "strand2",
            ]
        ]

    return df_sum.assign(score=lambda df: df["score"] / sum(total_counts))


def merge_pdf(
    pdf_generator: Iterable,
    merge_file: os.PathLike,
):
    pdf_files = []
    with pypdf.PdfWriter() as pdf_writer:
        for pdf_file in pdf_generator:
            pdf_writer.append(pdf_file)
            pdf_files.append(pdf_file)

        pdf_writer.write(merge_file)

    for pdf_file in pdf_files:
        pdf_file.unlink()


def draw_color_bar(
    cmap: str, vmin: float, vmax: float, label: str, outfile: os.PathLike
) -> None:
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)

    fig, ax = plt.subplots(figsize=(6, 1))
    fig.colorbar(mappable, cax=ax, orientation="horizontal", label=label)
    fig.tight_layout()
    fig.savefig(os.fspath(outfile))
    plt.close(fig)


class BlatSplice:
    def __init__(self, cfg: dict) -> None:
        self.databases = {
            assemble: f"{cfg[assemble]['2bit']}:{cfg[assemble]['chrom']}:{cfg[assemble]['start']}-{cfg[assemble]['end']}"
            for assemble in ["hg19", "mm10"]
        }
        self.blat = sh.Command("blat")

    def __call__(self, input: str, assemble: str) -> pd.DataFrame:
        # https://ucsc.crg.eu/FAQ/FAQblat.html for the parameter settings
        result = self.blat(
            "-out=blast8",
            "-stepSize=5",
            "-repMatch=2253",
            "-minScore=0",
            "-minIdentity=0",
            self.databases[assemble],
            "stdin",
            "stdout",
            _in=input,
        )

        return pd.read_csv(
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


def donor_acceptor_to_strand(donor: pd.Series, acceptor: pd.Series) -> pd.Series:
    return (
        pd
        .Series(data=["."] * len(donor))
        .where((donor != "GT") | (acceptor != "AG"), "+")
        .where((donor != "CT") | (acceptor != "AC"), "-")
    )


def read_interact(interact_file: os.PathLike) -> pd.DataFrame:
    return pd.read_csv(
        interact_file,
        sep="\t",
        names=[
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "value",
            "exp",
            "color",
            "sourceChrom",
            "sourceStart",
            "sourceEnd",
            "sourceName",
            "sourceStrand",
            "targetChrom",
            "targetStart",
            "targetEnd",
            "targetName",
            "targetStrand",
        ],
    )


def read_pairs(pairs_file: os.PathLike) -> pd.DataFrame:
    return pd.read_csv(
        pairs_file,
        sep="\t",
        skiprows=1,
        names=[
            "readID",
            "chrom1",
            "pos1",
            "chrom2",
            "pos2",
            "strand1",
            "strand2",
        ],
    )


def interact2pairs(interact_file: os.PathLike, pairs_file: os.PathLike) -> None:
    df = read_interact(interact_file)
    df = df.assign(
        pos1=lambda df: np.minimum(df["sourceEnd"], df["targetEnd"]) + 1,
        pos2=lambda df: np.maximum(df["sourceStart"], df["targetStart"]) + 1,
        strand1=lambda df: df["sourceStrand"].where(
            df["sourceStart"] <= df["targetStart"], df["targetStrand"]
        ),
        strand2=lambda df: df["targetStrand"].where(
            df["sourceStart"] <= df["targetStart"], df["sourceStrand"]
        ),
        chrom1=lambda df: df["sourceChrom"].where(
            df["sourceStart"] <= df["targetStart"], df["targetChrom"]
        ),
        chrom2=lambda df: df["targetChrom"].where(
            df["sourceStart"] <= df["targetStart"], df["sourceChrom"]
        ),
    ).rename(
        columns={
            "name": "readID",
        }
    )[
        [
            "readID",
            "chrom1",
            "pos1",
            "chrom2",
            "pos2",
            "strand1",
            "strand2",
        ]
    ]

    with open(pairs_file, "w") as fd:
        fd.write("## pairs format v1.0\n")
        df.to_csv(fd, sep="\t", header=False, index=False)


def pairs2bedpe(
    pairs_file: os.PathLike, bedpe_file: os.PathLike, total_count: int
) -> None:
    df = read_pairs(pairs_file)
    df = (
        df
        .rename(
            columns={
                "pos1": "end1",
                "pos2": "end2",
            }
        )
        .assign(
            start1=lambda df: df["end1"] - 1,
            start2=lambda df: df["end2"] - 1,
        )
        .groupby(
            [
                "chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "strand1",
                "strand2",
            ],
            as_index=False,
        )
        .agg(
            name=pd.NamedAgg("readID", lambda se: "|".join(se.tolist())),
            score=pd.NamedAgg("readID", "count"),
        )
        .assign(
            score=lambda df, total_count=total_count: (
                df["score"] / total_count * 1_000_000
            )
        )[
            [
                "chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "name",
                "score",
                "strand1",
                "strand2",
            ]
        ]
    )

    df.to_csv(
        bedpe_file,
        sep="\t",
        header=False,
        index=False,
    )

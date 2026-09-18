import matplotlib
import pandas as pd
import pypdf
from plotnine import (
    aes,
    element_text,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_fill_gradient,
    scale_y_discrete,
    theme,
)

from . import check
from .common import get_cpcdh_intron
from .utils import SelectTotalCount, blocks_string2tuple, clone2assemble, clone2treat

matplotlib.use("agg")


def splice(cfg: dict, assemble: str) -> None:
    df_intro = get_cpcdh_intron(
        cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
    ).query("name.str.lower().str.startswith('pcdha')")

    df = (
        (
            pd
            .read_feather(cfg["data_dir"] / "result" / "reads.feather")
            .assign(assemble=lambda df: df["clone"].map(clone2assemble))
            .query("assemble == @assemble")
            .reset_index(drop=True)
        )
        .assign(**{
            f"splice.{name}": lambda df, chrom=chrom, start=start, end=end: df[
                "ref_blocks"
            ].map(
                lambda ref_blocks, chrom=chrom, start=start, end=end: check.connect(
                    list(blocks_string2tuple(ref_blocks)), chrom, start, end, "+"
                )
            )
            for chrom, start, end, name in zip(
                df_intro["chrom"], df_intro["start"], df_intro["end"], df_intro["name"]
            )
        })
        .assign(**{
            f"precursor.{name}": lambda df, chrom=chrom, start=start: df[
                "ref_blocks"
            ].map(
                lambda ref_blocks, chrom=chrom, start=start - cfg["cover_threshold"], end=start + cfg["cover_threshold"]: (
                    check.cover(
                        list(blocks_string2tuple(ref_blocks)), chrom, start, end, "+"
                    )
                )
            )
            for chrom, start, name in zip(
                df_intro["chrom"], df_intro["start"], df_intro["name"]
            )
        })
    )

    df = (
        df
        .groupby(by=["exp", "protein", "clone", "rep", "query_name"], as_index=False)
        .agg(**{
            f"{opt}.{name}": pd.NamedAgg(column=f"{opt}.{name}", aggfunc="any")
            for name in df_intro["name"]
            for opt in ["splice", "precursor"]
        })
        .copy()
    )

    select_total_count = SelectTotalCount(cfg)
    df = (
        df
        .assign(treat=lambda df: df["clone"].map(clone2treat))
        .groupby(["exp", "protein", "treat"], as_index=False)
        .agg(**{
            f"{opt}.{name}": pd.NamedAgg(column=f"{opt}.{name}", aggfunc="sum")
            for name in df_intro["name"]
            for opt in ["splice", "precursor"]
        })
        .copy()
        .assign(
            total_count=lambda df: [
                select_total_count(exp, protein, treat)
                for exp, protein, treat in zip(df["exp"], df["protein"], df["treat"])
            ]
        )
    )

    df = (
        df
        .melt(
            id_vars=["exp", "protein", "treat", "total_count"],
            value_vars=[
                f"{opt}.{name}"
                for name in df_intro["name"]
                for opt in ["splice", "precursor"]
            ],
            var_name="opt_intron",
            value_name="count",
        )
        .assign(
            opt=lambda df: df["opt_intron"].str.split(".", expand=True)[0],
            intron=lambda df: df["opt_intron"].str.split(".", expand=True)[1],
        )
        .pivot_table(
            values="count",
            index=["exp", "protein", "treat", "total_count", "intron"],
            columns="opt",
        )
        .reset_index()
    )

    df = (
        pd
        .merge(
            df,
            df_intro[["chrom", "start", "end", "name"]],
            how="left",
            left_on="intron",
            right_on="name",
            validate="many_to_one",
        )
        .sort_values(by=["exp", "protein", "treat", "start"], ignore_index=True)
        .drop(columns="intron")
    )

    df.to_csv(cfg["data_dir"] / "result" / f"{assemble}_splice.csv", index=False)


def read_start_around_exon_start(cfg: dict, assemble: str) -> None:
    cpcdh_file = cfg["data_dir"] / "result" / f"{assemble}_cpcdh.csv"
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    tsses = df_cpcdh.query("name.str.lower().str.startswith('pcdha')")[
        ["start", "name"]
    ].reset_index(drop=True)

    df = (
        pd
        .read_feather(cfg["data_dir"] / "result" / "reads.feather")
        .assign(assemble=lambda df: df["clone"].map(clone2assemble))
        .query("assemble == @assemble")
        .assign(
            treat=lambda df: df["clone"].map(clone2treat),
            exp_protein_treat=lambda df: (
                df["exp"] + "_" + df["protein"] + "_" + df["treat"]
            ),
        )
        .assign(
            read_start=lambda df: df["ref_blocks"].map(
                lambda ref_blocks: check.start(list(blocks_string2tuple(ref_blocks)))
            )
        )[["exp_protein_treat", "read_start"]]
        .value_counts()
        .reset_index()
    )

    df = df.assign(**{
        exon: lambda df, tss=tss: df["read_start"] - tss
        for tss, exon in zip(tsses["start"], tsses["name"])
    })

    tss_extend = cfg["tss_extend"]
    df = (
        df
        .melt(
            id_vars=["exp_protein_treat", "read_start", "count"],
            value_vars=tsses["name"].tolist(),
            var_name="exon",
            value_name="relative",
        )
        .query("relative >= -@tss_extend and relative <= @tss_extend")
        .reset_index(drop=True)
        .assign(
            exon=lambda df: pd.Categorical(
                df["exon"],
                categories=tsses["name"].tolist(),
                ordered=True,
            )
        )
    )

    select_total_count = SelectTotalCount(cfg)
    pdf_files = []
    with pypdf.PdfWriter() as pdf_writer:
        for exp_protein_treat in df["exp_protein_treat"].unique():
            df_slice = (
                df
                .query("exp_protein_treat == @exp_protein_treat")
                .reset_index(drop=True)[["exon", "relative", "count"]]
                .set_index(["exon", "relative"])
                .reindex(
                    index=pd.MultiIndex.from_product(
                        [
                            tsses["name"],
                            list(range(-tss_extend, tss_extend + 1)),
                        ],
                        names=["exon", "relative"],
                    ),
                    fill_value=0,
                )
                .reset_index()
            )

            total_count = select_total_count(*exp_protein_treat.split("_"))
            df_slice = df_slice.assign(
                RPM=lambda df, total_count=total_count: (
                    df["count"] / total_count * 1_000_000
                ),
                RPM_round=lambda df: df["RPM"].round(2),
            )

            pdf_file = (
                cfg["data_dir"] / "result" / f"{assemble}_{exp_protein_treat}.pdf"
            )
            (
                ggplot(df_slice, mapping=aes(x="relative", y="exon"))
                + geom_tile(aes(fill="RPM"), color="#000000")
                + geom_text(aes(label="RPM_round"), size=6)
                + scale_fill_gradient(low="#FFFFFF", high="#FF0000")
                + scale_y_discrete(limits=tsses["name"].tolist()[::-1])
                + theme(
                    axis_text_x=element_text(angle=90, ma="right"), figure_size=(20, 20)
                )
                + labs(title=exp_protein_treat, x="position")
            ).save(pdf_file)

            pdf_writer.append(pdf_file)
            pdf_files.append(pdf_file)

        pdf_writer.write(
            cfg["data_dir"] / "result" / f"{assemble}_read_start_around_tss.pdf"
        )

    for pdf_file in pdf_files:
        pdf_file.unlink()


def get_pairs(cfg: dict, assemble: str) -> None:
    def splice_pair(ele: str) -> list:
        blocks = ele.split(";")
        pairs = []
        for i in range(len(blocks) - 1):
            block_chrom, block_start, block_end, block_strand = blocks[i].split(":")
            next_block_chrom, next_block_start, next_block_end, next_block_strand = (
                blocks[i + 1].split(":")
            )
            assert block_chrom == next_block_chrom, "inconsistent chrom"
            chrom1 = block_chrom
            chrom2 = next_block_chrom
            if block_start < next_block_start:
                pos1 = int(block_end) + 1
                strand1 = block_strand
                pos2 = int(next_block_start) + 1
                strand2 = next_block_strand
            else:
                pos1 = int(next_block_end) + 1
                strand1 = next_block_strand
                pos2 = int(block_start) + 1
                strand2 = block_strand

            pairs.append(f"{chrom1}:{pos1}:{strand1}:{chrom2}:{pos2}:{strand2}")

        return pairs

    result_file = cfg["data_dir"] / "result" / f"{assemble}_reads.feather"
    df = pd.read_feather(result_file).assign(
        treat=lambda df: get_treat(df),
        exp_protein_treat=lambda df: (
            df["exp"] + "_" + df["protein"] + "_" + df["treat"]
        ),
    )

    df = (
        df
        .query("not is_shadow and blocks.str.contains(';')")
        .reset_index(drop=True)[["exp_protein_treat", "query_name", "blocks"]]
        .assign(blocks=lambda df: df["blocks"].map(splice_pair))
        .explode("blocks", ignore_index=True)
    )
    df = pd.concat(
        [
            df[["exp_protein_treat", "query_name"]].rename(
                columns={"query_name": "readID"}
            ),
            df["blocks"]
            .str.split(":", expand=True)
            .rename(
                columns={
                    0: "chrom1",
                    1: "pos1",
                    2: "strand1",
                    3: "chrom2",
                    4: "pos2",
                    5: "strand2",
                }
            )[["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]],
        ],
        axis=1,
    )

    (cfg["data_dir"] / "result" / "hic" / "pairs").mkdir(exist_ok=True, parents=True)
    for exp_protein_treat in df["exp_protein_treat"].unique():
        with open(
            cfg["data_dir"] / "result" / "hic" / "pairs" / f"{exp_protein_treat}.pairs",
            "w",
        ) as fd:
            fd.write("## pairs format v1.0\n")
            df.query("exp_protein_treat == @exp_protein_treat")[
                ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
            ].to_csv(fd, sep="\t", header=False, index=False)


def get_interact(cfg: dict, assemble: str) -> None:
    def splice_interact(ele: str) -> list:
        blocks = ele.split(";")
        interact = []
        for i in range(len(blocks) - 1):
            block_chrom, block_start, block_end, block_strand = blocks[i].split(":")
            next_block_chrom, next_block_start, next_block_end, next_block_strand = (
                blocks[i + 1].split(":")
            )
            assert block_chrom == next_block_chrom, "chrom is not consistent"

            chrom = block_chrom
            chromStart = min(block_start, next_block_start)
            chromEnd = max(block_end, next_block_end)
            sourceChrom = block_chrom
            sourceStart = block_start
            sourceEnd = block_end
            sourceStrand = block_strand
            targetChrom = next_block_chrom
            targetStart = next_block_start
            targetEnd = next_block_end
            targetStrand = next_block_strand

            interact.append(
                f"{chrom}:{chromStart}:{chromEnd}:{sourceChrom}:{sourceStart}:{sourceEnd}:{sourceStrand}:{targetChrom}:{targetStart}:{targetEnd}:{targetStrand}"
            )

        return interact

    result_file = cfg["data_dir"] / "result" / f"{assemble}_reads.feather"
    df = pd.read_feather(result_file).assign(
        treat=lambda df: get_treat(df),
        exp_protein_treat=lambda df: (
            df["exp"] + "_" + df["protein"] + "_" + df["treat"]
        ),
    )

    df = (
        df
        .query("not is_shadow and blocks.str.contains(';')")
        .reset_index(drop=True)[["exp_protein_treat", "query_name", "blocks"]]
        .assign(blocks=lambda df: df["blocks"].map(splice_interact))
        .explode("blocks", ignore_index=True)
    )

    df = (
        pd
        .concat(
            [
                df[["exp_protein_treat", "query_name"]].rename(
                    columns={"exp_protein_treat": "exp", "query_name": "name"}
                ),
                df["blocks"]
                .str.split(":", expand=True)
                .rename(
                    columns={
                        0: "chrom",
                        1: "chromStart",
                        2: "chromEnd",
                        3: "sourceChrom",
                        4: "sourceStart",
                        5: "sourceEnd",
                        6: "sourceStrand",
                        7: "targetChrom",
                        8: "targetStart",
                        9: "targetEnd",
                        10: "targetStrand",
                    }
                ),
            ],
            axis=1,
        )
        .assign(score=0, color=0, value=1.0, sourceName=".", targetName=".")[
            [
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
            ]
        ]
        .sort_values(by=["chrom", "chromStart"], ignore_index=True)
    )

    (cfg["data_dir"] / "result" / "hic" / "interact").mkdir(exist_ok=True, parents=True)
    for exp_protein_treat in df["exp"].unique():
        df.query("exp == @exp_protein_treat").to_csv(
            cfg["data_dir"]
            / "result"
            / "hic"
            / "interact"
            / f"{exp_protein_treat}.bed",
            sep="\t",
            header=False,
            index=False,
        )

# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% vscode={"languageId": "shellscript"}
# !source src/check_splice/star.sh
# !index_ribosome
# !filter_ribosome_all
# !star_map_all


# %%
from check_splice import common, config, utils

cfg = config.pcdh()
utils.get_sample_bam(cfg).assign(
    total_count=lambda df: [common.get_bam_read_count(file) for file in df["file"]]
).to_csv(cfg["data_dir"] / "result" / "total_count.csv", index=False)

# %%
from check_splice import common, config, cpcdh

cfg = config.pcdh()
common.prepare_gene_bed12(
    gtffile=cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gtf",
    outfile=cfg["data_dir"] / "result" / "hg19.12.bed",
    addtional_filter=cpcdh.gene_bed12_hg19_cpcdh_filter,
)
common.prepare_gene_bed12(
    gtffile=cfg["data_dir"] / "data" / "mm10.ncbiRefSeq.gtf",
    outfile=cfg["data_dir"] / "result" / "mm10.12.bed",
    addtional_filter=cpcdh.gene_bed12_mm10_cpcdh_filter,
)

# %%
from check_splice import common, config, cpcdh

cfg = config.pcdh()
common.get_cpcdh_exon(
    gtffile=cfg["data_dir"] / "data" / "hg19.ncbiRefSeq.gtf", chrom="chr5"
).to_csv(cfg["data_dir"] / "result" / "hg19_cpcdh.csv", index=False)
cpcdh.fix_mm10_Pcdhgb8_CDS(
    common.get_cpcdh_exon(
        gtffile=cfg["data_dir"] / "data" / "mm10.ncbiRefSeq.gtf", chrom="chr18"
    )
).to_csv(cfg["data_dir"] / "result" / "mm10_cpcdh.csv", index=False)

# %%
from check_splice import config, sam

cfg = config.pcdh()
sam.merge_bam(cfg)

# %%
from check_splice import config, sam

cfg = config.pcdh()
sam.filter_bam(cfg)

# %%
from check_splice import config, sam

cfg = config.pcdh()
with open(cfg["data_dir"] / "result" / "reads.csv", "w") as fd:
    fd.writelines((f"{line}\n" for line in sam.parse_strand_sensitive_bam(cfg)))

# %%
from check_splice import config, sam

cfg = config.pcdh()
sam.group_read_blocks(cfg)

# %%
from check_splice import config, stat

cfg = config.pcdh()
stat.splice(cfg, "hg19")
stat.splice(cfg, "mm10")

# %%
from check_splice import config, draw

cfg = config.pcdh()
draw.splice_heatmap(cfg, "hg19")
draw.splice_heatmap(cfg, "mm10")

# %%
from check_splice import config, stat

cfg = config.pcdh()
stat.read_start_around_exon_start(cfg, "hg19")
stat.read_start_around_exon_start(cfg, "mm10")

# %%
from check_splice import config, stat

cfg = config.pcdh()
stat.get_exon_pre(cfg)

# %%
from check_splice import config, ply

cfg = config.pcdh()
ply.get_query_from_exp_protein_treat_query_name(cfg, "hg19")
ply.get_query_from_exp_protein_treat_query_name(cfg, "mm10")

# %%
from check_splice import config, ply

cfg = config.pcdh()
ply.get_plotly_interact(cfg, "hg19")
ply.get_plotly_interact(cfg, "mm10")

# %%
from check_splice import config, ply

cfg = config.pcdh()
ply.draw_interact(cfg, "hg19")
ply.draw_interact(cfg, "mm10")

# %%
from check_splice import config, interact

cfg = config.pcdh()
interact.get_interact(cfg)

# %%
from check_splice import config, interact

cfg = config.pcdh()
interact.interact_to_pairs(cfg)

# %%
from check_splice import config, interact

cfg = config.pcdh()
interact.pairs_to_bedpe(cfg)

# %%
from check_splice import config, interact

cfg = config.pcdh()
interact.diff_bedpe(cfg)

# %%
from check_splice import config, interact

cfg = config.pcdh()
interact.sum_bedpe(cfg)

# %%
from check_splice import config, interact

cfg = config.pcdh()
interact.filter_non_cpcdh_junction(cfg)

# %%
from check_splice import bigwig, config

cfg = config.pcdh()
bigwig.construct_artifact_bw(cfg, "hg19")
bigwig.construct_artifact_bw(cfg, "mm10")

# %%
from check_splice import bigwig, config

cfg = config.pcdh()
bigwig.construct_diff_bw(cfg)

# %%
from check_splice import cb, common, config

cfg = config.pcdh()
common.merge_pdf(
    cb.draw_links(cfg, cluster="alpha"),
    cfg["data_dir"] / "result" / "hic" / "draw" / "links.pdf",
)

# %%
from check_splice import cb, common, config

cfg = config.pcdh()
common.merge_pdf(
    cb.draw_pre_exons(cfg, cluster="alpha"),
    cfg["data_dir"] / "result" / "hic" / "draw" / "pre_exons.pdf",
)

# %%
from check_splice import cb, common, config

cfg = config.pcdh()
common.merge_pdf(
    cb.draw_reads(cfg),
    cfg["data_dir"] / "result" / "hic" / "draw" / "reads.pdf",
)

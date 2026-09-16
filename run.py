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
#     display_name: check-splice (3.13.11.final.0)
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
from check_splice import config, utils

cfg = config.pcdh()
utils.get_total_count(cfg)

# %%
from check_splice import config, cpcdh

cfg = config.pcdh()
cpcdh.prepare_gene_bed12(cfg)

# %%
from check_splice import config, cpcdh

cfg = config.pcdh()
cpcdh.get_cpcdh(cfg)

# %%
from check_splice import config, pcbs

cfg = config.pcdh()
pcbs.get_pCBS(cfg)

# %%
import json

from check_splice import config, process_all

cfg = config.pcdh()
for assemble in ["hg19", "mm10"]:
    with open(cfg["data_dir"] / "result" / f"{assemble}_reads.jsonl", "w") as fd:
        fd.writelines((f"{json.dumps(info)}\n" for info in process_all(cfg, assemble)))

# %%
from check_splice import config, utils

cfg = config.pcdh()
for assemble in ["hg19", "mm10"]:
    utils.jsonl2feather(
        jsonl_file=cfg["data_dir"] / "result" / f"{assemble}_reads.jsonl",
        feather_file=cfg["data_dir"] / "result" / f"{assemble}_reads.feather",
    )

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
stat.inrange_end_around_exon_end(cfg, "hg19")
stat.inrange_end_around_exon_end(cfg, "mm10")

# %%
from check_splice import config, sam

cfg = config.pcdh()
sam.merge_bam(cfg, "hg19")
sam.merge_bam(cfg, "mm10")

# %%
from check_splice import config, sam

cfg = config.pcdh()
sam.filter_bam_all(cfg)

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
from check_splice import config, stat

cfg = config.pcdh()
stat.get_pairs(cfg, "hg19")
stat.get_pairs(cfg, "mm10")

# %%
from check_splice import config, stat

cfg = config.pcdh()
stat.get_interact(cfg, "hg19")
stat.get_interact(cfg, "mm10")

# %%
from check_splice import cb, config

cfg = config.pcdh()
cb.pairs_to_bedpe(cfg)

# %%
from check_splice import cb, config

cfg = config.pcdh()
cb.diff_bedpe_all(cfg)

# %%
from check_splice import cb, config

cfg = config.pcdh()
cb.get_exon_pre(cfg, "hg19")
cb.get_exon_pre(cfg, "mm10")

# %%
from check_splice import cb, config

cfg = config.pcdh()
cb.construct_artifact_bw(cfg, "hg19")
cb.construct_artifact_bw(cfg, "mm10")

# %%
from check_splice import cb, config

cfg = config.pcdh()
cb.construct_diff_bw(cfg, "hg19")
cb.construct_diff_bw(cfg, "mm10")

# %%
from check_splice import cb, config

cfg = config.pcdh()
cb.draw_all(cfg, "hg19")
cb.draw_all(cfg, "mm10")

# %%
from check_splice import cb, config

cfg = config.pcdh()
cb.draw_reads_all(cfg, "hg19")
cb.draw_reads_all(cfg, "mm10")

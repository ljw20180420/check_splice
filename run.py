#!/usr/bin/env python

from check_splice import config, utils

cfg = config.pcdh()
utils.jsonl2feather(
    jsonl_file=cfg["data_dir"] / "result" / "reads.jsonl",
    feather_file=cfg["data_dir"] / "result" / "reads.feather",
)

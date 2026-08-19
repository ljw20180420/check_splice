#!/usr/bin/env python

import pandas as pd

from check_splice import config

cfg = config.pcdh()
df = pd.read_json(cfg["data_dir"] / "result" / "reads.jsonl", lines=True)
df.to_feather(cfg["data_dir"] / "result" / "reads.feather")

#!/usr/bin/env python

import json

from check_splice import config, process_all

cfg = config.pcdh()
with open(cfg["data_dir"] / "result" / "reads.jsonl", "w") as fd:
    fd.writelines((f"{json.dumps(info)}\n" for info in process_all(cfg)))

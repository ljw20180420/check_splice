#!/usr/bin/env python

from check_splice import config, stat

cfg = config.pcdh()
stat.inrange_end_around_exon_end(cfg)

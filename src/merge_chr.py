#!/usr/bin/env python3

import os
import sys
import gzip
from mod_report import HEADER

header_len = len(HEADER.split('\n'))-1

fnames = sys.argv[1:]; print(fnames)
outfn = os.path.join(os.path.dirname(fnames[0]), "mod.gz")
with gzip.open(outfn, "wt") as out:
    for fi, fn in enumerate(fnames):
        for li, line in enumerate(gzip.open(fn, "rt")):
            if fi and li<header_len:
                continue
            out.write(line)

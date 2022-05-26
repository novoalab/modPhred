#!/usr/bin/env python3
# Combine intermedaite mod.gz.*.gz files. 
# USAGE: src/modPhred/src/merge_chr.py outdir/mod.gz.chr*.gz

import os
import sys
import gzip
from mod_report import HEADER

header_len = len(HEADER.split('\n'))-1

fnames = sys.argv[1:]
print("Merging mod.gz files for %s chromosomes..."%len(fnames))
outfn = os.path.join(os.path.dirname(fnames[0]), "mod.gz")
with gzip.open(outfn, "wt") as out:
    i = 0
    for fi, fn in enumerate(fnames):
        for li, line in enumerate(gzip.open(fn, "rt")):
            if fi and li<header_len:
                continue
            out.write(line)
            i += 1
print(" %s lines merged to %s!"%(i, outfn))

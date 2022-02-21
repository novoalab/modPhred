#!/usr/bin/env python3

import os, sys, resource
import numpy as np
from datetime import datetime

VERSION = '1.1a'

# set max 3 modifications per each base - this is done automatically depending on the model
# but is minimally three
#MaxModsPerBase = 3
MaxProb = 256
QUALS = np.array(list(map(chr, range(33, 127)))) #"".join(map(chr, range(33, 127)))
# it's only DNA as in SAM U should be A
base2complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

HEADER = """###
# Welcome to modPhred (ver. %s)!
# 
# Executed with: %s
#
# For each bam file 4 values are stored for every position:
# - depth of coverage (only positions with >=%s X in at least one sample are reported)
# - accuracy of basecalling (fraction of reads having same base as reference, ignoring indels)
# - frequency of modification (fraction of reads with modification above given threshold)
# - median modification probability of modified bases (0-1 scaled). 
#
# If you have any questions, suggestions or want to report bugs,
# please use https://github.com/novoalab/modPhred/issues.
# 
# Let's begin the fun-time with Nanopore modifications...
###
chr\tpos\tref_base\tstrand\tmod\t%s
"""

def make_counter(count=0):
    """Simple counter"""
    def inner():
        nonlocal count
        count += 1
        return count
    return inner

def warning(info, m=3, counter=make_counter()):
    """Return at most m of warning messages"""
    # update counter
    count = counter()#; print(count)
    # skip if too many warnings returned
    if m and count>m: return
    if info: logger(info, add_timestamp=0, add_memory=0)
    if count==m:
        logger("Future warnings will be skipped!", add_timestamp=0, add_memory=0)

def memory_usage(childrenmem=True, div=1024.):
    """Return memory usage in MB including children processes"""
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / div
    if childrenmem:
        mem += resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / div
    return mem
        
def logger(info, add_timestamp=1, add_memory=1, out=sys.stderr):
    """Report nicely formatted stream to stderr"""
    info = info.rstrip('\n')
    memory = timestamp = ""
    if add_timestamp:
        timestamp = "[%s]"%str(datetime.now()).split(".")[0] #"[%s]"%datetime.ctime(datetime.now())
    if add_memory:
        memory = " [mem: %5.0f MB]"%memory_usage()
    out.write("%s %s%s\n"%(timestamp, info, memory))


def get_alphabet(output_alphabet, mods, canonical_bases="ACGTU", force_rna=0):
    """Return correctly ordered alphabet of bases with modified bases
    and 2 dictionaries: symbol2modbase & canonical2mods.
    """
    if isinstance(mods, str): mods = mods.split() # list of base_mod_long_names ie 6mA 5mC
    else: mods = mods[:]
    output_alphabet = list(output_alphabet) # AYCZGT
    # get ordered alphabet and dictionaries
    symbol2modbase, canonical2mods = {}, {}
    alphabet, _mods = [], []
    # decide if DNA or RNA # ABCGH4TW
    if force_rna: output_alphabet[output_alphabet.index("T")] = "U"
    bases = "".join(b for b in output_alphabet if b in canonical_bases)
    canonical2mods = {b: [] for b in bases}
    idx = -1 
    for b in output_alphabet:
        alphabet.append(b)
        if b in bases:
            idx += 1
        else:
            if b.isdigit():
                canonical2mods[bases[idx+1]].append(b)
            else:
                canonical2mods[bases[idx]].append(b)
            symbol2modbase[b] = mods.pop(0)
    # get base2positions
    base2positions = {}
    idx = 0
    for b in bases:
        base2positions[b] = [idx, ]
        idx += 1
        for mb in canonical2mods[b]:
            base2positions[b].append(idx)
            idx += 1
    #print(bases, alphabet, symbol2modbase, canonical2mods, base2positions)
    return alphabet, symbol2modbase, canonical2mods, base2positions

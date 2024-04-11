#!/usr/bin/env python3
desc="""Plot density of modification probability.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Torredembarra, 21/12/2023
"""

import os, pysam, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from datetime import datetime

def get_strand2positionFn(outdir, bed, strands="+-"):
    strand2positionFn = {s: os.path.join(outdir, "pos.{}.txt".format(s)) for s in strands}
    strand2out = {s: open(strand2positionFn[s], "wt") for s in strands}
    for line in open(bed, 'rt'):
        chrom, s, e, name, score, strand = line[:-1].split('\t')[:6]
        if strand in strand2out:
            strand2out[strand].write("{} {}\n".format(chrom, e))
    for a, out in strand2out.items():
        out.close()
    return strand2positionFn

def get_quals_for_regions(bams, strand2positionFn, boi, mapq=10, baseq=0):
    args=["--ff", "3840", "--no-output-ins", "--no-output-ins", "--no-output-del", "--no-output-del", 
          "--no-output-ends", "-q%s"%mapq, "-Q%s"%baseq, "-B", *bams]
    
    bam2quals = [[] for bam in bams]
    for strand, fn in strand2positionFn.items():
        mpileup = pysam.mpileup("-l", fn, *args)
        boi = boi.upper() if strand=="+" else boi.lower()
        for line in mpileup.rstrip().split('\n'):
            ldata = line.split('\t')
            if len(ldata) < 3:
                continue
            for bami, (bases, quals) in enumerate(zip(ldata[4::3], ldata[5::3])):
                sel = np.array(list(bases))==boi
                bam2quals[bami].append("".join(np.array(list(quals))[sel]))
    return bam2quals

def plot(outdir, mer, probs, N_positions, labels=("A", "m6A"), vlines=[], ext="pdf"):
    fig1, ax1 = plt.subplots(figsize=(7, 5))
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    for label, _probs in zip(labels, probs):
        ax1.hist(_probs, bins=31, range=(0, 1), label=label, alpha=0.5)
        sns.kdeplot(_probs, label=label, ax=ax2)#, bw_adjust=0.1)
    ax1.set_ylabel("Number of reads")
    ax2.set_ylabel("Density")
    for ax in (ax1, ax2):
        ax.set_xlabel("Modification probability")
        ax.set_title("{} : {:,} reads from {:,} positions".format(mer, len(probs[0]), N_positions))
        ax.legend()
    # trim X-axis to 0-1
    ax2.set_xlim(0, 1)
    if vlines:
        ax1.vlines(vlines, *ax1.get_ylim(), ls=":", color="black", alpha=0.33)
        ax2.vlines(vlines, *ax2.get_ylim(), ls=":", color="black", alpha=0.33)
    fig1.savefig(os.path.join(outdir, "{}.hist.{}".format(mer, ext)))
    fig2.savefig(os.path.join(outdir, "{}.density.{}".format(mer, ext)))

def bam2density(outdir, bams, bed, mer, samples=[], max_reads=100,
                vlines=[], log=sys.stderr):
    
    log.write("Creating output directory...\n")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if not samples:
        samples = [os.path.basename(bam) for bam in bams]
    log.write("{} BAM file(s) for samples: {}\n".format(len(bams), ", ".join(samples)))
    # prepare 
    log.write("Parsing BED file...\n")
    strand2positionFn = get_strand2positionFn(outdir, bed)
    # 
    boi = mer[len(mer)//2]
    log.write("Processing BAM file(s)...\n")
    quals = get_quals_for_regions(bams, strand2positionFn, boi)
    N_positions = len(quals[0])
    quals_merged = ["".join(q[:max_reads] for q in quals[i]) for i in range(len(quals))]
    # get probs: quals are 0-30 scaled
    probs = [(np.array(list(quals_merged[i])).view(np.int32) - 33) / 30 for i in range(len(quals))]
    log.write("Plotting...\n")
    plot(outdir, mer, probs, N_positions, samples, vlines)

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-o", "--outdir", required=True, help="output directory")
    parser.add_argument("-b", "--bed", required=True, help="regions BED")
    parser.add_argument("-i", "--input", nargs="+", help="input BAM file(s)")
    parser.add_argument("-s", "--samples", nargs="+", help="sample name(s)")
    parser.add_argument("-m", "--mer", default="DRACH",
                        help="mer of interest [%(default)s]. Central base will be used as base of interest")
    parser.add_argument("-q", "--mapq", default=15, type=int,
                        help="min mapping quality [%(default)s]")
    parser.add_argument("-d", "--maxDepth", default=100, type=int,
                        help="max depth of coverage [%(default)s]")
    parser.add_argument("--vlines", default=[], type=float, nargs="+", 
                        help="max depth of coverage [%(default)s]")
    o = parser.parse_args()

    bam2density(o.outdir, o.input, o.bed, o.mer, o.samples, o.maxDepth, o.vlines)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)

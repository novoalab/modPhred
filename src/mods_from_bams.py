#!/usr/bin/env python3
desc="""Report positions likely being modified from mapped reads
with modification probabilities encoded as base qualities. 

The pipeline is performed in following steps:
- mod_encode.py: encodes modifications probabilities from Fast5 as base qualities in FastQ
- mod_align.py: aligns reads from FastQ files to reference and stores BAM files for each sample
- mod_report.py: reports modified sites from BAM files that fulfil given filtering criteria
- mod_plot.py: performs QC and plotting

More info at: https://github.com/lpryszcz/modPhred
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 3/02/2023
"""

import os, sys
from datetime import datetime
from mod_report import load_info, mod_report, logger, VERSION


def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument("-i", "--bams", nargs="+", help="input directories with Fast5 files")
    parser.add_argument("-o", "--outdir", default="modPhred", help="output directory [%(default)s]")
    parser.add_argument("-f", "--fasta", required=1, type=argparse.FileType('r'), help="reference FASTA file")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="min mapping quality [%(default)s]")
    parser.add_argument("-d", "--minDepth", default=25, type=int, help="min depth of coverage [%(default)s]")
    parser.add_argument("--minModFreq", default=0.05, type=float, help="min modification frequency per position [%(default)s]")
    parser.add_argument("--minModProb", default=0.50, type=float, help="min modification probability per base [%(default)s]")
    parser.add_argument("-b", "--bed", help="BED file with regions to analyse [optionally]")
    parser.add_argument("--chr", default=[], nargs="+", help="chromosome(s) to process [all]")
    parser.add_argument("--step", type=int, default=20000, help="region size [%(default)s]")
    parser.add_argument("-t", "--threads", default=6, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-u", "--unspliced", action='store_true', help="don't use spliced alignment for RNA")
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    # we need name, not file object here
    o.fasta = o.fasta.name
    
    logger("===== Welcome, welcome to modPhred pipeline! =====", add_memory=0)
    if os.path.isdir(o.outdir):
        logger("Output directory exits. Steps completed previously will be skipped!")
    
    # process everything only if outfile does not exists
    outfn = os.path.join(o.outdir, "mod.gz")
    if o.chr:
        outfn += "."+",".join(o.chr)+".gz"
    if not os.path.isfile(outfn):
        # load info
        data = load_info(o.outdir)
        MaxPhredProb = data["MaxPhredProb"]
        bamfiles = o.bams #data["bam"]
        #bamfiles.sort()
        # index
        logger("Indexing bam file(s)...")
        for fn in bamfiles:
            if not os.path.isfile(fn+".bai"):
                cmd = "samtools index %s"%fn
                if o.verbose:
                    sys.stderr.write(" %s\n"%cmd)
                os.system(cmd)
        # BAM > modifications
        # get can2mods ie {'A': ['6mA'], 'C': ['5mC'], 'G': [], 'T': []}
        can2mods = {b: [data["symbol2modbase"][m] for m in mods]
                    for b, mods in data["canonical2mods"].items()}
        # U>T patch for RNA mods
        if "U" in can2mods:
            can2mods["T"] = can2mods["U"]
        mod_report(outfn, bamfiles, o.fasta, o.threads, o.bed, MaxPhredProb, can2mods, 
                   o.mapq, o.minDepth, o.minModFreq, o.minModProb,
                   o.chr, o.step, logger=logger)
    else:
        logger(" %s exists!\n"%outfn)

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

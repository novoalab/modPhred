#!/usr/bin/env python3
desc="""Get regions for mer. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Torredembarra, 21/12/2023
"""

import os, pysam, re, sys
from datetime import datetime

base2complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

def get_regions_for_mer(outfn, fasta, mers, strands="+-", score="."):
    bed = "{}\t{}\t{}\t{}\t{}\t{}\n"
    with open(outfn, "wt") as out:
        faidx = pysam.FastaFile(fasta)
        for ref in faidx.references:
            seq = faidx.fetch(ref)
            for mer in mers:
                offset = len(mer)//2
                # forward matches
                if "+" in strands:
                    for m in re.compile(mer).finditer(seq):
                        s = m.start()+offset
                        out.write(bed.format(ref, s, s+1, mer, score, "+"))
                # reverse complement matches
                if "-" in strands:
                    merrc = "".join(base2complement[b] for b in mer[::-1])
                    for m in re.compile(merrc).finditer(seq):
                        s = m.start()+offset
                        out.write(bed.format(ref, s, s+1, mer, score, "-"))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-o", "--outfn", required=True, help="output file name")
    parser.add_argument("-f", "--fasta", required=True, help="FastA file")
    parser.add_argument("-m", "--mer", nargs="+", default=["GGACT", ],
                        help="mer of interest [%(default)s]. Central base will be used as base of interest")
    parser.add_argument("-s", "--strands", default="+-",
                        help="strand(s) to report [%(default)s]")
    o = parser.parse_args()

    get_regions_for_mer(o.outfn, o.fasta, o.mer, o.strands)

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



#!/usr/bin/env python3
desc="""Add base qualities to SAM entries
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 29/06/2021
"""

import os, sys, pysam
from datetime import datetime
from guppy_encode import VERSION

def get_rname_qual(fastq):
    """Return """
    for fq in fastq:
        for r in pysam.FastqFile(fq):
            yield r.name, r.quality #r.get_quality_array()

def sam_add_qual(instream, outstream, fastq, tag, logger=sys.stderr.write):
    """Add original base qualities to SAM entries"""
    quals = get_rname_qual(fastq)
    sam = pysam.AlignmentFile(instream)
    out = pysam.AlignmentFile(outstream, "w", header=sam.header) # wb is ~2x slower!
    pname = ""
    for i, a in enumerate(sam, 1):
        # get next qual for each individual read
        if a.qname!=pname:
            rname, qual = next(quals)
            pname = a.qname
            if rname != pname:
                logger("[WARNING] Read name mismatch @ %s alignment: %s vs %s\n"%(i, rname, pname))
        if not a.is_unmapped:
            # reverse qual if alignment is reverse
            if a.is_reverse:
                qual = qual[::-1]
            # clip the qual if alignment is hard-clipped
            s = a.cigar[0][1] if a.cigar[0][0]==5 else 0
            e = -a.cigar[-1][1] if a.cigar[-1][0]==5 else None
            qual = qual[s:e]
        # add to alignment & write
        a.set_tag(tag, qual)
        out.write(a)
    out.close()

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument("-i", "--input", default=sys.stdin, help="input stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, help="output stream [stdout]")
    parser.add_argument("-f", "--fastq", nargs="+", help="input FastQ files")
    parser.add_argument("-t", "--tag", default="OQ", help="SAM tag to store qualities [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # get alignments
    sam_add_qual(o.input, o.output, o.fastq, o.tag)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    #sys.stderr.write("#Time elapsed: %s\n"%dt)

#!/usr/bin/env python3
desc="""Use mapping with shared buffer to save BAM files from multiple inputs

This helps to keep low memory footprint for concurent mapping of 
multiple multi_fast5 files, especially for large reference like human. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 22/07/2021
"""

import gzip, sys, mappy, pysam
from datetime import datetime
from multiprocessing import Pool

def get_alignments_from_mappy(fastq, header, MD=True, cs=False):
    """Generator of alignments from mappy. 

    Reference index has to be already loaded!
    """
    # get index buffer from shared memory
    thr_buf = mappy.ThreadBuffer()
    # process reads
    for name, seq, qual in mappy.fastx_read(fastq):
        # iterate hits
        primary_count = 0
        for hit in aligner.map(seq, buf=thr_buf, MD=MD, cs=cs): 
            # catch secondary alignments
            if hit.is_primary:
                primary_count += 1
            a = mappy2sam(name, seq, qual, hit, header, primary_count)
            yield a

def init_args(*args):
    """Share globals with pool of workers"""
    global aligner
    aligner, = args

def mappy2sam(name, seq, qual, h, header, primary_count=0):
    """Create AlignedSegment and add it to sam file"""
    a = pysam.AlignedSegment(header)
    skip_st = "%sS"%h.q_st if h.q_st else ""
    end_skip = len(seq)-h.q_en
    skip_en = "%sS"%end_skip if end_skip else ""
    if h.strand<0: # or h.trans_strand<0
        a.is_reverse = True
        seq = mappy.revcomp(seq)
        qual = qual[::-1]
        skip_st, skip_en = skip_en, skip_st
    a.qname = name
    a.reference_name = h.ctg
    a.reference_start = h.r_st
    a.seq = seq
    a.qual = qual
    a.mapq = h.mapq
    a.cigarstring = skip_st + h.cigar_str + skip_en
    a.is_secondary = False if h.is_primary else True
    a.is_supplementary = True if h.is_primary and primary_count>1 else False
    a.set_tags([("NM", h.NM, "i"), ("MD", h.MD, "Z")]) # ("AS", ,"i"),
    return a
    
def worker(args):
    """Worker performing actual mapping"""
    fastq, ref = args
    # get FastA index
    faidx = pysam.FastaFile(ref)#; faidx.references, faidx.lengths
    # get output bam
    bam = "%s.bam"%fastq
    sam = pysam.AlignmentFile(bam, mode="wb", reference_names=faidx.references, 
                          reference_lengths=faidx.lengths)
    # get index buffer from shared memory
    thr_buf = mappy.ThreadBuffer()
    # process reads
    for name, seq, qual in mappy.fastx_read(fastq):
        # iterate hits
        primary_count = 0
        for hit in aligner.map(seq, buf=thr_buf, MD=True): #, cs=False
            # catch secondary alignments
            if hit.is_primary:
                primary_count += 1
            a = mappy2sam(name, seq, qual, hit, sam.header, primary_count)
            sam.write(a)
    sam.close()
    return bam

def mapping(fasta, fastq, rna, threads, log=sys.stderr):
    """Main mapping function"""
    # init aligner
    log.write("Loading index from %s...\n"%fasta)
    aligner = mappy.Aligner(fasta, preset="spliced" if rna else "map-ont", k=13)
    args = [(f, fasta) for f in fastq]
    # start pool of workers initialising shared aligner
    log.write("Initialising pool of %s workers...\n"%threads)
    p = Pool(threads, initializer=init_args, initargs=(aligner, ))
    # align individual FastQ files
    log.write("Aligning %s files...\n"%len(fastq))
    for i, bam in enumerate(p.imap_unordered(worker, args), 1):
        log.write(" %s %s  \r"%(i, bam))
    p.close()
    
# -t3 
# 0:00:49 56,980 map-ont
# 0:00:51 57,391 map-ont k=13

# -t1
# 0:02:06 57,391 map-ont k=13

# 0:03:37 -t3 HUMAN 10G memory, 20G at peak during index creation
# 0:03:25 minimap2 --MD -ax spliced -t 3 -k 13 $ref  | samtools view -Sb > test.bam
# 0:16:33 time for f in _archives/raw/xPore/HEK293T-WT-rep1/batch_?.fast5.bam.fq.gz; do minimap2 --MD -ax spliced -t 3 -k 13 $ref $f | samtools view -Sb > $f.bam2; done

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    #parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("-i", "--fastq", nargs="+", help="FastA/Q file(s)")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-t", "--threads", default=3, type=int, help="number of cores to use (and connections to accept) [%(default)s]")
    parser.add_argument("-l", "--log", default=sys.stderr, type=argparse.FileType("w"), help="logging destination [stderr]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    mapping(o.fasta, o.fastq, o.rna, o.threads, o.log)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s    \n"%dt)

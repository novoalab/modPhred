#!/usr/bin/env python3
desc="""Report positions likely being modified from mapped reads
with modification probabilities encoded as base qualities. 

The pipeline is performed in following steps:
- mod_encode.py: basells Fast5, align and encode modifications probabilities as base qualities in BAM
- mod_report.py: reports modified sites from BAM files that fulfil given filtering criteria
- mod_plot.py: performs QC and plotting

Dependencies: numpy mappy pandas pysam seaborn pyguppyclient (matching guppy version)

More info at: https://github.com/lpryszcz/modPhred
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 21/06/2019
"""

import glob, gzip, os, pysam, sys, zlib
from datetime import datetime
from multiprocessing import Pool
import numpy as np
import pandas as pd
from guppy_encode import *
from mod_plot import load_bed
from pathlib import Path
import mod_plot

alphabet = "NACGT" # first base should be always N!!!
base2index = {b: i for i, b in enumerate(alphabet)}
for i, b in enumerate(alphabet.lower()):
    base2index[b] = i

# CIGAR operations
"""Op BAM Description +1Q +1R
M 0 alignment match (can be a sequence match or mismatch) yes yes
I 1 insertion to the reference yes no
D 2 deletion from the reference no yes
N 3 skipped region from the reference no yes
S 4 soft clipping (clipped sequences present in SEQ) yes no
H 5 hard clipping (clipped sequences NOT present in SEQ) no no
P 6 padding (silent deletion from padded reference) no no
= 7 sequence match yes yes
X 8 sequence mismatch yes yes
    """
def _match(refi, readi, bases): return refi+bases, readi+bases, True
def _insertion(refi, readi, bases): return refi, readi+bases, False
def _deletion(refi, readi, bases): return refi+bases, readi, False
def _skip(refi, readi, bases): return refi, readi, False
code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                 2: _deletion, 3: _deletion, 4: _insertion, 5: _skip}

def baseq2float(b, q, minModProb, MaxPhredProb, base2index=base2index):
    """Return float representation of base and it's modification probability

    For example A (base 1) with modifications probability 31 will be stored as 1.31.
    Ignore probabilities below minModPorb for speed. 
    """
    # ignore mod prob below minModProb
    if (q%MaxPhredProb)<minModProb*MaxPhredProb:
        return base2index[b]
    # q is max 93 so it should be padded with 0 for q<10
    # ie q=1 should be 0.01 and q=10 would be 0.10
    return float('{}.{:02d}'.format(base2index[b], q))

def store_blocks(a, start, end, i, quals, minModProb, MaxPhredProb):
    """Store base calls from aligned blocks. """
    readi, refi = 0, a.pos
    for ci, (code, bases) in enumerate(a.cigar):
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        # skip if current before start
        if refi<=start:
            continue
        # typical alignment part
        if data:
            if prefi<start:
                bases -= start-prefi
                preadi += start-prefi
                prefi = start
            if refi>end:
                bases -= refi-end
            if bases<1:
                break
            for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases])):
                #if b in base2index: 
                # add qual - 28b per int, 50b per 1 char and 24b for float
                quals[i][prefi-start+ii].append(baseq2float(b, q, minModProb, MaxPhredProb))
    return quals

def is_qcfail(a, mapq=15, flag=1792): #
    """Return True if alignment fails QC

    By default skips:
    - alignment with mapping quality below 15
    - not primary alignment (256)
    - read fails platform/vendor quality checks (512)
    - read is PCR or optical duplicate (1024)
    - and alignments without base qualities
    
    Bit flags are in brackets from https://broadinstitute.github.io/picard/explain-flags.html.
    """
    if a.mapq<mapq or a.flag&flag or not a.query_qualities:
        return True
    return False

def bam2calls(sam, ref, start, end, mapq, minModProb, MaxPhredProb):
    """Generator of basecalls and mod qualities from BAM file encoded as floats for +/- strand"""
    # prepare storage for quals -- this can be memory & speed optimised
    quals = [[[] for x in range(end-start+1)], [[] for x in range(end-start+1)]]
    p = 0
    empty = []
    for a in sam.fetch(ref, start, end):
        # skip low quality alignments, not primary, QC fails, duplicates or supplementary algs
        if is_qcfail(a, mapq): continue
        # yield some results
        while p<a.pos-start:
            yield quals[0][p]
            yield quals[1][p]
            # and release a bit of memory
            quals[0][p] = quals[1][p] = empty
            p += 1
        # get transcript strand: this is for secondstrand or dRNAseq
        si = 1 if a.is_reverse else 0 # si==0: +/for; si==1: -/rev
        # store alignment blocks
        quals = store_blocks(a, start, end, si, quals, minModProb, MaxPhredProb)
    # yield last bit of calls
    while p<len(quals[0]):
        yield quals[0][p]
        yield quals[1][p]
        p += 1
    
def fasta2bases(faidx, ref, start, end):
    """Generator of individual bases from FastA file.

    The output consists of: 
    - position in reference (1-based)
    - strand integer (0 for plus, 1 for minus)
    - strand as +/-
    - base (complement for -)
    - base index in alphabet (relative to + strand)
    """
    if ref not in faidx:
        raise StopIteration
    seq = faidx.fetch(ref, start, end).upper()
    for pos, refbase in enumerate(seq, start+1):
        # update refbase and get appropriate cols
        if refbase in base2index:
            refi = base2index[refbase]
        else:
            sys.stderr.write("[WARNING] %s:%s %s not in alphabet (%s). Marking as N.\n"%(ref, pos, refbase, alphabet))
            refbase, refi = "N", 0 # N
        yield pos, 0, "+", refbase, refi
        yield pos, 1, "-", base2complement[refbase], refi

def get_modCount(strand_calls, bams, maxNmods, MaxPhredProb, minModProb):
    """Return modCount array that count for each bam and every base
    number of modified and unmodified bases"""
    modCount = np.zeros((len(bams), len(alphabet), maxNmods+1), dtype='uint64')
    probs = [[[[] for mi in range(maxNmods)] for base in alphabet] for bam in bams]
    for bami, _calls in enumerate(strand_calls):
        for bf in _calls:
            # store as unmodified if modification prob below minModQual (reported as int)
            if isinstance(bf, int):
                bi, mi = bf, -1
            else:
                # get base as int
                bi = int(bf)
                # get modification as phred quality, which mod it is and store its probability
                mq = (100.*(bf-int(bf)))
                mi = int(mq/MaxPhredProb) # which mod is it (first, second, third etc. last here is unmodified)
                mp = mq%MaxPhredProb / (MaxPhredProb-1)
                # store modprob
                probs[bami][bi][mi].append(mp)
            modCount[bami, bi, mi] += 1
    return modCount, probs
        
def get_calls(position, fasta, bams, MaxPhredProb, can2mods,
              mapq=15, minDepth=25, minModFreq=0.01, minModProb=0.5):
    """Return modified positions from given region"""
    # get region info
    ref, start, end = position
    info = []
    posinfo = "%s\t%s\t%s>%s%s"
    # get actual max no of modifications per base
    maxNmods = max(map(len, can2mods.values()))
    # init faidx and sam once and reuse - this saves ~0.5s per region
    ## to it via init_args and global args
    #faidx = pysam.FastaFile(fasta)
    #sams = [pysam.AlignmentFile(bam) for bam in bams]
    refparser = fasta2bases(faidx, ref, start, end)
    parsers = [bam2calls(sam, ref, start, end, mapq, minModProb, MaxPhredProb) for sam in sams]
    for data in zip(refparser, *parsers):
        (pos, si, strand, refbase, refi) = data[0]
        calls = data[1:]
        modCount, probs = get_modCount(calls, bams, maxNmods, MaxPhredProb, minModProb)
        # skip if low coverage in all samples
        coverage = modCount.sum(axis=2).sum(axis=1)
        if coverage.max() < minDepth: continue
        baseacc = 1.* modCount.sum(axis=2)[:,refi] / coverage # simply number of As if ref A divided by total coverage for given positions
        freq = modCount[:,:,:-1] / coverage[:, None, None]
        # report all - ignore base 0 (N)
        for bi, base in enumerate(alphabet[1:], 1):
            if si: base = base2complement[base]
            if not can2mods[base]: continue
            # iterate mods for given base
            for mi in range(maxNmods):
                # skip if low freq of modification at given position
                if np.nanmax(freq[:, bi, mi])<minModFreq: continue
                #print(bi, base, np.max(freq[:,bi]), freq[:,bi,mi]); print(probs[0][bi][mi], probs[0][bi][-1])
                # get normalised mean mod quality taking into account that there may be more than 1 mod type
                medqual = map(np.median, (probs[fi][bi][mi] for fi in range(len(bams))))
                mod = can2mods[base][mi] # mod always is relevant to the read base, not ref
                text = "\t".join([ref, str(pos), refbase, strand, mod, 
                                  "\t".join("%s\t%.3f\t%.3f\t%.3f"%d for d in zip(coverage, baseacc, freq[:,bi, mi], medqual))])
                #print(text)
                info.append(text + "\n")
    return "".join(info)

def worker(args):
    # ignore all warnings
    import warnings
    warnings.filterwarnings('ignore')
    data = get_calls(*args)
    return data
    
def get_consecutive(data, stepsize=1):
    """Return consecutive windows allowing given max. step size"""
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)

def get_coverage_for_ref(sam, ref, mapq, reflen): 
    """Return coverage from sam (pysam.AlignmentFile)"""
    coverage = np.zeros(reflen, dtype='uint16')
    for a in sam.fetch(ref):
        if is_qcfail(a, mapq): continue
        if a.blocks:
            for s, e in a.blocks: coverage[s:e] += 1
        else: coverage[a.pos:a.aend] += 1
    return coverage

def get_covered_regions(bams, fasta, threads=6, mapq=20, minCov=10, maxdist=16000, step=20000):
    """Return regions covered by at least minCov reads.
    
    Here, we only analyse the file with lower number of reads.
    """
    # iterate chromosomes
    sams = [pysam.AlignmentFile(bam, threads=threads) for bam in bams]
    # file with fewer alignemnts first
    sams = list(sorted(sams, key=lambda x: x.mapped))
    sam = sams[0]
    faidx = pysam.FastaFile(fasta)
    regions = []
    logger("Retrieving regions covered by %s+ reads..."%minCov)
    # first run idxstats and skip references with <10 reads altogether - usefull for transcript alignements
    ref2algs = {s.contig: s.mapped for s in sam.get_index_statistics()}
    ref2len = {r: l for r, l in zip(sam.references, sam.lengths)}
    bases = 0
    # get references with enough alignments
    refs = [r for r, c in ref2algs.items() if c>=minCov and r in faidx]
    # catch transcript alignments
    if len(ref2len)>3000:
        # just get transcripts with minCov alignments
        for ref in refs:
            s, e = 0, ref2len[ref]
            bases += e
            # further split regions for max windows
            while s < e-step:
                regions.append((ref, s, s+step))
                s += step
            regions.append((ref, s, e))
    else:
        for ri, ref in enumerate(refs, 1):
            sys.stderr.write(" %s / %s %s ...\r"%(ri, len(sam.references), ref))
            #if ri>100:break
            # get min coverage from bam file that has fewer mapped reads
            coverage = get_coverage_for_ref(sam, ref, mapq, ref2len[ref])
            # get regions with minCov coverage
            covered = np.where(coverage>=minCov)[0]
            for positions in get_consecutive(covered, maxdist):
                if len(positions)<1: continue
                s, e = positions[0]+1, positions[-1]+1
                bases += e-s
                # further split regions for max windows
                while s < e-step:
                    regions.append((ref, s, s+step))
                    s += step
                regions.append((ref, s, e))
    logger(" {:,} bases in {:,} regions to process.\n".format(bases, len(regions)))
    return regions

def init_args(*args):
    """Share globals with pool of workers"""
    global sams, faidx
    bams, fasta = args
    sams = [pysam.AlignmentFile(bam) for bam in bams]
    faidx = pysam.FastaFile(fasta)

def mod_report(outfn, bam, fasta, threads, regionsfn, MaxPhredProb, can2mods,
               mapq=15, minDepth=25, minModFreq=0.1, minModProb=0.5, logger=sys.stderr.write):
    """Get modifications from bam files"""
    logger("Reporting positions that are likely modified to %s ..."%outfn)
    # write header
    out = gzip.open(outfn, "wt")
    cols = ["depth", "basecall_accuracy", "mod_frequency", "median_mod_prob"]
    out.write(HEADER%(VERSION, " ".join(sys.argv), minDepth,
                      "\t".join(" ".join((b, c)) for b in bam for c in cols)))
    # get covered regions - only for RNA!!!
    logger(" Getting regions covered by at least %s reads..."%minDepth)
    if regionsfn:
        regions = load_bed(regionsfn)
    else:
        regions = get_covered_regions(bam, fasta, threads, mapq, minDepth)
    logger("  %s regions to process..."%len(regions))
    # define imap, either pool of processes or map
    p = Pool(threads, initializer=init_args, initargs=(bam, fasta), maxtasksperchild=100)
    # start genotyping
    i = 0
    parser = p.imap(worker, ((pos, fasta, bam, MaxPhredProb, can2mods,
                              mapq, minDepth, minModFreq, minModProb)
                             for pos in regions))
    for i, data in enumerate(parser, 1):
        sys.stderr.write(" %s / %s [memory: %7.1f Mb]\r"%(i, len(regions), memory_usage()))
        out.write(data)
    # finalise
    #logger(" %s regions processed"%i)
    out.close()

def get_rgb(mod_color, freq):
    """Return R,G,B color definition"""
    rgb = mod_color * int(round(255*freq))
    return ",".join(map(str, rgb))

def mod_bed(data, outdir, symbol2modbase, minModFreq=0.1):
    """Generate bedMethyl (BED-formatted) modification output:
    
    Reference chromosome or scaffold
    Start position in chromosome
    End position in chromosome
    Name of item (short modification name)
    Score (median mod prob) from 0-1000. 
    Strandedness, plus (+), minus (-), or unknown (.)
    Start of where display should be thick (start codon)
    End of where display should be thick (stop codon)
    Color value (RGB) - color strengh depends on modification frequency; different color to various modifications, although if more than 7 mods the colors may repeat
    Coverage, or number of reads
    Percentage of reads that show given modification at this position in the genome/transcriptome
    """
    # RGB colors: use G, B, R
    colors = np.array([(0, 1, 0), (0, 0, 1), (1, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)], dtype=int)
    # load info & get mod2color - colors may repeat
    mod2color = {}
    for i, m in enumerate(symbol2modbase.values()):
        mod2color[m] = colors[i%len(colors)]
        
    outfn = os.path.join(outdir, "mod.bed")
    logger("Saving modified positions with max frequency as %s (bedMethyl file) ..."%outfn)
    if os.path.isfile(outfn):
        logger(" file %s exists!\n"%outfn)
        return
    # replace NaN with 0
    data.replace(to_replace=np.nan, value=0, inplace=True) #data.fillna(0, inplace=True) - this replaces entire rows
    # save bed for every input BAM file and 1 combined
    out = open(outfn, "w") 
    # get metrics to save score
    depthcols = list(filter(lambda x: x.endswith("depth"), data.columns))
    probcols = list(filter(lambda x: x.endswith("median_mod_prob"), data.columns))
    freqcols = list(filter(lambda x: x.endswith("mod_frequency"), data.columns))
    # get output files
    fnames = [x[:-len("mod_frequency")-1] for x in freqcols]#; print(fnames)
    outs = [open(fn+".bed", "w") for fn in fnames]
    # info ans save
    logger(" and separately for every BAM file as %s/*.bed ..."%os.path.dirname(fnames[0]))
    columns = ["chr", "pos", "ref_base", "strand", "mod"] + depthcols + probcols + freqcols
    line = "%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%i\n"
    for i, r in data[columns].iterrows():
        # save combined output with max values for score, depth & freq
        out.write(line%(r.chr, r.pos-1, r.pos, r["mod"], round(1000*r[probcols].max()), r.strand,
                        r.pos-1, r.pos,  get_rgb(mod2color[r["mod"]],r[freqcols].max()), 
                        r[depthcols].max(), round(100*r[freqcols].max()))) 
        # and scores for each individual sample
        for f, depth, score, freq in zip(outs, r[depthcols], r[probcols], r[freqcols]):
            if freq<minModFreq: continue
            f.write(line%(r.chr, r.pos-1, r.pos, r["mod"], round(1000*score), r.strand,
                          r.pos-1, r.pos, get_rgb(mod2color[r["mod"]], freq), depth, round(100*freq)))
    out.close()
    
def get_can2mods(md):
    """Return can2mods"""
    alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(md['base_mod_alphabet'], md['base_mod_long_names'])
    # get can2mods ie {'A': ['6mA'], 'C': ['5mC'], 'G': [], 'T': []}
    can2mods = {b: [symbol2modbase[m] for m in mods]
                for b, mods in canonical2mods.items()}
    # U>T patch for RNA mods
    if "U" in can2mods: can2mods["T"] = can2mods["U"]
    return can2mods
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument("-i", "--indirs", nargs="+", help="input directories with Fast5 files")
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("--sensitive", action='store_true', help="use sensitive mapping parameters ie tRNA")
    parser.add_argument("-o", "--outdir", default="modPhred", help="output directory [%(default)s]")
    parser.add_argument("-f", "--fasta", required=1, type=argparse.FileType('r'), help="reference FASTA file")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="min mapping quality [%(default)s]")
    parser.add_argument("-d", "--minDepth", default=25, type=int, help="min depth of coverage [%(default)s]")
    parser.add_argument("--minModFreq", default=0.05, type=float, help="min modification frequency per position [%(default)s]")
    parser.add_argument("--minModProb", default=0.50, type=float, help="min modification probability per base [%(default)s]")
    parser.add_argument("-b", "--bed", help="BED file with regions to analyse [optionally]")
    parser.add_argument("-t", "--threads", default=6, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("--tag", default="", help="SAM tag to store original qualities ie. OQ [skipped]")
    guppy = parser.add_argument_group("Basecalling options")
    guppy.add_argument("-c", "--config", default="dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg", help="guppy model [%(default)s]")
    guppy.add_argument("--host", "--guppy_basecall_server", default=None,
                        help="guppy server hostname or path to guppy_basecall_server binary [assumes files are already basecalled with modifications]")
    guppy.add_argument("--port", default=5555, type=int,
                        help="guppy server port (this is ignored if binary is provided) [%(default)s]")
    guppy.add_argument("--device", default="cuda:0", help="CUDA device to use (works only if --host guppy_basecall_server_path) [%(default)s]")
    guppy.add_argument("--timeout", default=20*60, type=int, help="timeout in seconds to process each Fast5 file [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    # we need name, not file object here
    o.fasta = o.fasta.name
    
    logger("===== Welcome, welcome to modPhred pipeline! =====", add_memory=0)
    if os.path.isdir(o.outdir):
        logger("Output directory exits. Steps completed previously will be skipped!")
    
    # encode modifications in BAM
    bamfiles = mod_encode(o.outdir, o.indirs, o.fasta, o.threads, o.rna, o.sensitive, 
                          o.config, o.host, o.port, o.recursive, o.device, o.timeout, o.tag)
        
    # load info
    fnames, basecount, mods2count, md = get_mod_data(bamfiles[0])
    alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(md['base_mod_alphabet'], md['base_mod_long_names'])
    MaxPhredProb = md["MaxPhredProb"]
    can2mods = get_can2mods(md)
    
    # process everything only if outfile does not exists
    outfn = os.path.join(o.outdir, "mod.gz")
    if not os.path.isfile(outfn):
        # BAM > modifications
        # the index is already created by mod_encode, but let's keep it
        logger("Indexing BAM file(s)...")
        for fn in bamfiles:
            if not os.path.isfile(fn+".bai") and not os.path.isfile(fn+".csi"):
                logger(" %s"%fn)
                pysam.index(fn)
        mod_report(outfn, bamfiles, o.fasta, o.threads, o.bed, MaxPhredProb, can2mods, 
                   o.mapq, o.minDepth, o.minModFreq, o.minModProb, logger=logger)
    else:
        logger(" %s exists!\n"%outfn)

    # load data
    logger("Loading modification data...")
    data = pd.read_csv(outfn, sep="\t", header=len(HEADER.split('\n'))-2, index_col=False)
    # plot
    logger("Plotting...")
    mod_plot.mod_plot(outfn, data=data)
    # get bed-formatted output
    mod_bed(data, o.outdir, symbol2modbase, o.minModFreq)
    # plot pairwise plots - those will affect column names therefore has to be done at the end
    try:
        # skip pairplots if more than 4 samples
        if len(o.indirs)>4:
            logger("More than 4 samples: pairwise plots skipped. To generate them, execute:")
            logger(" src/mod_plot.py --scatter -i %s"%outfn)
        elif len(o.indirs)>1:
            mod_plot.plot_scatter(outfn, data=data)
        # nothint to plot if single sample
    except Exception as e:
        logger("[ERROR] Plotting scatterplots failed (%s).\n"%str(e))

    logger("All finished! Have a nice day :)")

if __name__=='__main__': 
    t0 = datetime.now()
    os.setpgrp() # create new process group, become its leader    
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except Exception as err:
        import signal, traceback
        sys.stderr.write(traceback.format_exc()+"\n")
        os.killpg(0, signal.SIGTERM) # terminate all processes in my group
    finally:
        dt = datetime.now()-t0
        sys.stderr.write("#Time elapsed: %s    \n"%dt)


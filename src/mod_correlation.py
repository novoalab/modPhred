#!/usr/bin/env python3
desc="""Report correlation between all modified positions

More info at: https://github.com/lpryszcz/modPhred
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 3/02/2020
"""

import glob, gzip, os, pysam, sys, zlib
from collections import Counter
from datetime import datetime
from multiprocessing import Pool
from  matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from guppy_encode import *
from mod_report import is_qcfail, base2index, code2function

complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

def baseq2int(b, q, minModProb, MaxPhredProb, base2index=base2index, maxMods=10):
    """Return int representation of modification

    For example:
    - A (base 1) with modification probability of
      - 31 will be stored as 1 (first mod of A) because 1 + 0*maxMods
    - C (base 2) with modification probability of
      - 50 will be stored as 12 (second mod of C) because 2 + 1*maxMods
    """
    # ignore mod prob below minModProb - store as 255
    if (q%MaxPhredProb)<minModProb*MaxPhredProb:
        return 255
    # story only mod per base
    return 1 + q//MaxPhredProb + maxMods*(base2index[b]-1)

def store_blocks(a, start, end, pos2idx, quals, readidx, minModProb, MaxPhredProb):
    """Store base calls from aligned blocks. """
    readi, refi = 0, a.pos
    stored = False
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
            # process bases of aligned block
            for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases]), prefi+1):
                if ii in pos2idx:
                    # get complement if read mapped to rev strand
                    if a.is_reverse: b = complement[b]
                    # store most likely modification as int if larger than minModProb
                    modi = baseq2int(b, q, minModProb, MaxPhredProb)
                    if modi:
                        #if modi<255: print(readidx, ii, pos2idx[ii], b, q, base2index[b], modi)
                        #if b=="C" and 1000<ii<1050: print(readidx, ii, pos2idx[ii], b, q, base2index[b], modi)
                        quals[pos2idx[ii]][readidx] = modi
                        stored = True
    if stored:
        readidx += 1
    return quals, readidx

def bam2calls(bam, ref, positions, mapq, minModProb, MaxPhredProb, maxDepth=25000):
    """Generator of basecalls and mod qualities from BAM file encoded as floats for +/- strand"""
    sam = pysam.AlignmentFile(bam)
    start, end = positions[0]-1, positions[-1]+1
    # prepare storage
    pos2idx = {p: i for i, p in enumerate(positions)}
    quals = np.zeros((len(positions), maxDepth), dtype='uint8') # up to 64 mods per base!
    # stop if ref not in sam file
    if ref not in sam.references:
        if ref.startswith('chr') and ref[3:] in sam.references:
            ref = ref[3:]
        elif 'chr%s'%ref in sam.references:
            ref = 'chr%s'%ref
    posidx = readidx = 0#; print(bam)
    for a in sam.fetch(ref, start, end):
        # skip low quality alignments, not primary, QC fails, duplicates or supplementary algs
        if is_qcfail(a, mapq):
            continue
        # yield some results and release a bit of memory
        while positions[posidx]<a.pos:
            # yield only reads informative for given position
            yield quals[:, quals[posidx]>0] #np.all((quals[posidx]>0, quals[posidx]<255), axis=0)]
            posidx += 1
        # store modifications from alignment blocks
        quals, readidx = store_blocks(a, start, end, pos2idx, quals, readidx, minModProb, MaxPhredProb)
        # release some more memory
        if readidx==maxDepth:
            #sys.stderr.write("[INFO] maxDepth reached\n"); break ### to be done
            first = np.argwhere(quals[posidx, :readidx])[0][0]
            if first<0.05*maxDepth: first = int(0.05*maxDepth)
            quals = np.hstack((quals[:, first:], np.zeros((len(positions), first), dtype='uint8')))
            sys.stderr.write("[INFO][%s] Dropped %s previous columns/reads\n"%(bam, first))
            readidx -= first
    # yield last bit of calls
    while posidx<len(positions):
        # yield only reads informative for given position
        yield quals[:, quals[posidx]>0]
        posidx += 1

def chr2modcorr(outfn, bams, region, chrdata, mapq, mindepth, minModProb, MaxPhredProb, minmodreads=10):
    """Calculate correlation between modifications"""
    cols = ["chr", "pos", "mod", "strand"]
    # get chr and positions
    ref, positions = chrdata.chr.unique()[0], np.unique(chrdata.pos.to_numpy())
    corrs = np.zeros((len(positions), len(positions)), dtype="float32")
    corrs[:] = np.nan
    logger(" %s with %s modified positions > %s"%(region, len(positions), outfn))
    parsers = [bam2calls(bam, ref, positions, mapq, minModProb, MaxPhredProb) for bam in bams]
    for i, calls in enumerate(zip(*parsers)):
        # stack all reads - those are already prefiltered for only those modified for given position
        calls = np.hstack(calls)
        sys.stderr.write(" %s  \r"%i)
        # get modified positions
        #mod = calls!=255
        # get positions with mindepth
        enoughdepth = np.where(np.sum(calls>0, axis=1)>=mindepth)[0]
        print(Counter(calls[i]), enoughdepth.sum())
        # store correlations between this positions
        for j in filter(lambda x: x>=i, enoughdepth):
            # mod in i and j and take balanced number of modified reads for each position
            modi = np.argwhere(np.all((calls[i]>0, calls[i]<255, calls[j]>0), axis=0)) #calls[i]>0, calls[i]<255 #calls[i]==1
            modj = np.argwhere(np.all((calls[j]>0, calls[j]<255, calls[i]>0), axis=0)) #calls[j]>0, calls[j]<255 #calls[j]==1
            lessmod = min(len(modi), len(modj))
            sel = np.unique(list(modi[:lessmod]) + list(modj[:lessmod]))
            if len(sel)<minmodreads: continue
            '''# take only those reads that are modified in any of the samples
            sel = np.all((np.any(mod[[i, j]], axis=0), calls[j]>0), axis=0)
            # skip if less than 10 reads with modification at least in 1 position
            if sel.sum()<minmodreads:
                continue'''
            # get those that are modified in both - unmod are 255
            same = calls[i, sel] == calls[j, sel]
            corr = 2*(np.mean(same)-0.5) if np.any(same) else -1 #; print(chrdata[chrdata.pos==positions[j]][cols].to_numpy(), i, j, corr, same.sum(), sel.sum())
            corrs[i, j] = corrs[j, i] = corr        
        '''# get positions modified for all positions  
        modreads = np.all((calls>0, calls<255), axis=0)
        modreadsum = modreads.sum(axis=1) 
        # get positions with mindepth
        enoughdepth = np.where(np.sum(calls>0, axis=1)>=mindepth)[0]
        # store correlations between this positions
        for j in filter(lambda x: x>=i, enoughdepth):
            # choose _i that it always has more modified positions than _j
            if np.argmax(modreadsum[[i, j]]):
                _j, _i = i, j
            else:
                _i, _j = i, j
            # select only reads modified in the position with more modifications
            # and with bases called in both
            sel = np.all((modreads[_i], calls[_j]>0), axis=0) #; print(i, j, modreadsum[[i, j]])
            # skip if less than 10 reads with modification at least in 1 position
            if sel.sum()<minmodreads:
                continue
            # get those that are modified in both - unmod are 255
            same = calls[_j, sel]!=255 # calls[i, sel] == calls[j, sel]
            corr = 2*(np.mean(same)-0.5) if np.any(same) else -1 #; print(chrdata[chrdata.pos==positions[j]][cols].to_numpy(), i, j, corr, same.sum(), sel.sum())
            corrs[i, j] = corrs[j, i] = corr'''
        #return
    #print(corrs[:10,:10]); return
    # store array and
    np.savetxt(outfn, corrs, fmt='%.3f', header=",".join(map(str, positions)), delimiter=',',
               footer=",".join("%s%s"%(m, s) for m, s in zip(chrdata["mod"], chrdata.strand)))
    return corrs

def get_data_for_regions(data, regions):
    """Return list of DataFrames limited to given regions"""
    regionsData = []
    for region in regions:
        chrom, s, e = region, 0, 0
        # get strand chr:start-end
        strand = None
        if region.endswith(("+", "-")):
            region, strand = region[:-1], region[-1]
        # get chr, start & end from chr:start-end
        if ":" in region and "-" in region:
            chrom, se = region.split(':')
            s, e = map(int, se.split("-"))
        df = data[data.chr==chrom]
        if e:
            df = df[np.all((df.pos>=s, df.pos<=e), axis=0)]
        if strand:
            df = df[df.strand==strand]
        regionsData.append(df)
    return regionsData
        
def mod_correlation(infn, ext="png", logger=logger, data=False, overwrite=False, regions=[], samples=[], 
                    minfreq=0.20, mindepth=10, minModProb=0.5, mapq=15, strand=None, mod=None):
    outdir = os.path.join(os.path.dirname(infn), "correlations")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # get bamfiles
    bamfiles = glob.glob(os.path.join(os.path.dirname(infn), "minimap2", "*.bam"))
    bamfiles.sort()
    # load info
    fnames, basecount, mods2count, md = get_mod_data(bamfiles[0])
    alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(md['base_mod_alphabet'], md['base_mod_long_names'])
    MaxPhredProb = md["MaxPhredProb"]
    # U>T patch for RNA mods
    if "U" in can2mods:
        can2mods["T"] = can2mods["U"]
    
    # parse data
    if isinstance(data, bool):
        logger("Loading %s ...\n"%infn)
        data = pd.read_csv(infn, sep="\t", header=len(HEADER.split('\n'))-2, index_col=False,
                           dtype={"chr": object, "pos": int}) # ADD TO ALL
    # filter by min freq and depth
    mfreqcols = list(filter(lambda x: x.endswith('mod_frequency'), data.columns)); mfreqcols
    depthcols = list(filter(lambda x: x.endswith('depth'), data.columns)); depthcols
    filters = [data.loc[:, mfreqcols].max(axis=1)>minfreq, data.loc[:, depthcols].max(axis=1)>mindepth]
    # add filters for strand and modification
    if mod:
        filters.append(data["mod"]==mod)
    data = data[np.all(filters, axis=0)]
    #print(data.shape, data.head())
    # limit by region AND CONSIDER LIMITING COV TO 2-3x median?
    if regions:
        # get regions logger(" limiting to %s regions: %s\n"%(len(regions), ",".join(regions)))
        regionsData = get_data_for_regions(data, regions)
        logger("Processing %s region(s): %s ...\n"%(len(regions), ",".join(regions)[:3]))
    else:
        # get chromosomes
        regions = data.chr.unique()
        #if strand: filters.append(data.strand==strand)
        regionsData = (data[data.chr==ref] for ref in regions)
        logger("Processing %s chromosome(s): %s ...\n"%(len(regions), ",".join(regions)[:3]))
    if data.shape[0]<1:
        logger("[mod_plot][ERROR]  %s row(s) found in %s\n"%(data.shape[0], infn))
        return
    # process regions/chromosomes
    for ref, chrdata in zip(regions, regionsData):
        # define output
        fn = "%s.csv.gz"%ref
        if mod: fn = "%s.%s.csv.gz"%(ref, mod)
        outfn = os.path.join(outdir, fn)
        if overwrite or not os.path.isfile(outfn):
            # generate data
            corrs = chr2modcorr(outfn, bamfiles, ref, chrdata, mapq, mindepth, minModProb, MaxPhredProb)
        else:
            # load data
            corrs = np.loadtxt(outfn, delimiter=",")
        # plot
        plot_heatmap(corrs, chrdata, ref, outfn, ext=ext)
        # maybe plot below freq? that would be cool, right?

def collapse_axes(xlab, ylab):
    """Return X/Y labels for unique X positions"""
    _xlab, _ylab, pos = [xlab[0]], [ylab[0]], [0]
    for i in range(1, len(xlab)):
        if xlab[i] == xlab[i-1]:
            _ylab[-1] += "; %s"%ylab[i]
        else:
            _xlab.append(xlab[i])
            _ylab.append(ylab[i])
            pos.append(i)
    return _xlab, _ylab, pos
        
def plot_heatmap(corrs, chrdata, ref, outfn, figsize=(12, 10), dim=100000, ext="svg", simpleY=True):
    """Plot heatmaps"""
    # narrow by strands
    xlab = chrdata.pos[:dim].to_numpy()
    ylab = ["%s %s"%(m, s) for m, s in zip(chrdata["mod"][:dim], chrdata.strand[:dim])]
    mod2count = Counter(chrdata["mod"])
    xlab, ylab, pos = collapse_axes(xlab, ylab)
    # use unique names on Y
    if simpleY:
        _ylab = [ylab[i] if ylab[i]!=ylab[i-1] else "" for i in range(1, len(ylab))]
        _ylab.insert(0, ylab[0])
        ylab = _ylab
    # switch axes labels
    #xlab, ylab = ylab, xlab
    logger(" Plotting %s modified positions in %s:%s-%s"%(len(chrdata), ref, xlab[0], xlab[-1]))
    #mask = np.zeros_like(corrs)
    #mask[np.triu_indices_from(mask)] = True
    #f, ax = plt.subplots(figsize=figsize)
    fig = plt.figure(figsize=figsize)
    # add title
    mcounts = "; ".join("%s: %s"%(m, c) for m, c in mod2count.items())
    fig.suptitle("\n%s\n%s modifications: %s"%(ref, len(chrdata), mcounts))
    fig.subplots_adjust(top=0.75)    
    with sns.axes_style("white"):
        ax = sns.heatmap(corrs, vmin=-1, vmax=1, center=0, cmap="RdBu_r", #mask=mask, #fmt="d"
                         xticklabels=xlab, yticklabels=ylab)
        #ax.xaxis.tick_top() 
    ax.set_xlabel("Modified positions")
    ax.set_ylabel("Modifications at those positions [with +/- strand]")
    # calculate global frequency
    depthcols = list(filter(lambda c: c.endswith("depth"), chrdata.columns))
    mfreqcols = list(filter(lambda x: x.endswith('mod_frequency'), chrdata.columns))
    mcountcols = ["%s modcount"%c.split()[0] for c in depthcols]
    for cc, fc, dc in zip(mcountcols, depthcols, mfreqcols):
        chrdata[cc] = chrdata[fc] * chrdata[dc]
    chrdata["avgfreq"] = chrdata[mcountcols].sum(axis=1) / chrdata[depthcols].sum(axis=1)
    ax2 = fig.add_axes([.125, 0.77, .62, .10], anchor="N", sharex=ax)
    ax2.bar(np.arange(len(corrs))+0.5, chrdata["avgfreq"].to_numpy()[pos])
    ax2.set_ylim((0, 1)); ax2.set_ylabel("mod\nfreq")
    ax2.xaxis.tick_top() #ax2.set_xticklabels([]) #
    ax2.set_xticklabels(xlab, rotation=90)
    fig.savefig(outfn+".%s"%ext) #show()
            
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--input", default="modPhred/mod.gz", help="input file [%(default)s]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="min mapping quality [%(default)s]")
    parser.add_argument("-d", "--minDepth", default=25, type=int, help="min depth of coverage [%(default)s]")
    parser.add_argument("--minModFreq", default=0.20, type=float, help="min modification frequency per position [%(default)s]")
    parser.add_argument("--minModProb", default=0.50, type=float, help="min modification probability per base [%(default)s]")
    #parser.add_argument("-s", "--strand", default=None, choices=["+", "-"], help="select strand [include both]")
    parser.add_argument("--mod", default="", help="filter only 1 modification [analyse all]")
    parser.add_argument("-r", "--regions", nargs="+", default=[], help="regions to process [all chromosomes]")    
    parser.add_argument("-w", "--overwrite", action="store_true", help="overwrite existing output")    
    parser.add_argument("-e", "--ext", default="svg", help="figure format/extension [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    if not o.regions:
        logger("Processing entire chromosomes - consider narrowing to certain regions!")
        
    mod_correlation(o.input, ext=o.ext, mapq=o.mapq, overwrite=o.overwrite, regions=o.regions, mod=o.mod, 
                    minfreq=o.minModFreq, mindepth=o.minDepth, minModProb=o.minModProb)
    logger("Finished\n")

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

#!/usr/bin/env python3
desc="""Report correlation between all modified positions

More info at: https://github.com/lpryszcz/modPhred
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 3/21/2020
"""

import glob, gzip, os, pickle, pysam, sys, zlib
from collections import Counter
from datetime import datetime
from multiprocessing import Pool
from  matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from guppy_encode import HEADER, VERSION, logger, memory_usage, load_info, base2complement, MaxModsPerBase
from mod_report import is_qcfail, base2index, code2function
from sklearn.decomposition import PCA
#from scipy.cluster.hierarchy import dendrogram, linkage

complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

def store_blocks(a, start, end, df, reads, minModProb, MaxPhredProb, can2mods, pos2base):
    """Store base calls from aligned blocks. """
    readi, refi, stored = 0, a.pos, 0
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
            # process bases of aligned block - here ii is 1-off because of prefi+1
            for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases]), prefi+1):
                if ii in pos2base:
                    strand = "+"
                    if a.is_reverse:
                        b = complement[b]
                        strand = "-"
                    if b!=pos2base[ii]:
                        continue
                    # get which modification it is and its probability
                    modi = q//MaxPhredProb
                    if (ii, can2mods[b][modi], strand) not in df.columns: continue
                    modprob = (q%MaxPhredProb)/(MaxPhredProb-1)#; print(ii, b, q, modi, modprob)
                    #df.loc[[len(reads)], "%s %s %s"%(ii, can2mods[b][modi], strand)] = modprob
                    df.loc[[len(reads)], (ii, can2mods[b][modi], strand)] = modprob
                    stored += 1
    return df, stored

def bam2calls(bams, chrdata, mapq, minModProb, MaxPhredProb, can2mods,
              maxDepth=100, minStored=0, minAlgFrac=0.3, strand=None):
    """Generator of basecalls and mod qualities from BAM file encoded as floats for +/- strand"""
    ref = chrdata.chr.iloc[0]
    start, end = chrdata.pos.iloc[0]-1, chrdata.pos.iloc[-1]+1
    pos2base = {r.pos: r.ref_base for idx, r in chrdata.iterrows() if can2mods[r.ref_base]}
    # prepare data frame with reads in rows x modified positions in columns
    columns = pd.MultiIndex.from_frame(chrdata[['pos', 'mod', "strand"]])
    df = pd.DataFrame(np.zeros((maxDepth*len(bams), chrdata.shape[0])), columns=columns)#; print(df.head())
    # add sample as last column
    df["sample"] = df["read_strand"] = ""
    logger(" %s:%s-%s %s with %s modified positions"%(ref, start, end, strand, len(pos2base)))
    # process bam files
    reads = []
    for bi, bam in enumerate(bams, 1):
        # get sample name modPhred/curlcakes/minimap2/RNA010220191_m5C.bam -> RNA010220191_m5C
        sample = os.path.basename(bam)[:-4]#; print(sample)
        sys.stderr.write("  %s / %s %s > %s              \r"%(bi, len(bams), bam, sample))
        readi = 0
        sam = pysam.AlignmentFile(bam)
        for a in sam.fetch(ref, start, end):
            # skip low quality alignments, not primary, QC fails, duplicates or supplementary algs
            # or algs shorter than 30% of the region
            s, e = max(start, a.pos), min(end, a.aend)
            if is_qcfail(a, mapq) or (e-s) < minAlgFrac*(end-start):
                continue
            # skip reads from wrong strand
            if strand=="+" and a.is_reverse or strand=="-" and not a.is_reverse:
                continue
            # store modifications from alignment blocks
            df, stored = store_blocks(a, start, end, df, reads, minModProb, MaxPhredProb, can2mods, pos2base)
            # keep read if enough mods
            if stored>=minStored:
                df.loc[len(reads), "read_strand"] = "-" if a.is_reverse else "+"
                reads.append(a.qname) #"%s: %s"%(sample, a.qname))
                readi += 1
            else:
                df.loc[len(reads)] = 0
            # process up to maxDepth of reads per BAM file
            if readi==maxDepth:
                break
        # store samples
        df.loc[len(reads)-readi:len(reads), "sample"] = sample #.split("_")[-1]
    # strip rows that have no reads & rename rows
    df = df.iloc[:len(reads)]
    df.index = reads
    return df

def get_data_for_regions(data, regions):
    """Return list of DataFrames limited to given regions"""
    if len(regions)==1 and os.path.isfile(regions[0]) or os.path.islink(regions[0]):
        fn = regions[0]
        regions = []
        for l in open(fn):
            ldata = l[:-1].split("\t")
            r = "%s:%s-%s"%tuple(ldata[:3])
            if len(ldata)>5 and ldata[5] in "+-":
                r += ldata[5]
            regions.append(r)
    #print(regions)
    regionsData, strands = [], []
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
        strands.append(strand)
    return regions, regionsData, strands
        
def mod_cluster(outdir, infn, bamfiles, args, data=False, ext="png", logger=logger, read_strand_lut = {"+": "r", "-": "b"}):
    """Cluster reads base on their modification profiles"""
    regions = args.regions
    if not outdir:
        outdir = os.path.join(os.path.dirname(infn), "clusters")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # load info
    moddata = load_info(os.path.dirname(infn))
    MaxPhredProb = moddata["MaxPhredProb"]
    if not bamfiles:
        bamfiles = moddata["bam"]
        bamfiles.sort()
    # BAM > modifications
    # get can2mods ie {'A': ['6mA'], 'C': ['5mC'], 'G': [], 'T': []}
    can2mods = {b: [moddata["symbol2modbase"][m] for m in mods]
                for b, mods in moddata["canonical2mods"].items()}
    if "U" in can2mods:
        can2mods["T"] = can2mods["U"]#; print(can2mods)
    #print(MaxPhredProb, can2mods, bamfiles)
    
    # parse data
    if isinstance(data, bool):
        logger("Loading %s ...\n"%infn)
        data = pd.read_csv(infn, sep="\t", header=len(HEADER.split('\n'))-2, index_col=False,
                           dtype={"chr": object, "pos": int}) # ADD TO ALL
    # filter by min freq and depth
    mfreqcols = list(filter(lambda x: x.endswith('mod_frequency'), data.columns)); mfreqcols
    depthcols = list(filter(lambda x: x.endswith('depth'), data.columns)); depthcols
    filters = [data.loc[:, mfreqcols].max(axis=1)>args.minfreq,
               data.loc[:, depthcols].max(axis=1)>args.mindepth]
    # add filters for strand and modification
    if args.mod:
        filters.append(data["mod"]==args.mod)
    data = data[np.all(filters, axis=0)]
    #print(data.shape, data.head())
    # limit by region AND CONSIDER LIMITING COV TO 2-3x median?
    if regions:
        # get regions logger(" limiting to %s regions: %s\n"%(len(regions), ",".join(regions)))
        regions, regionsData, strands = get_data_for_regions(data, regions)
        logger("Processing %s region(s): %s ...\n"%(len(regions), ",".join(regions)[:3]))
    else:
        # get chromosomes
        regions = data.chr.unique()
        strands = [None] * len(regions)
        #if strand: filters.append(data.strand==strand)
        regionsData = (data[data.chr==ref] for ref in regions)
        logger("Processing %s chromosome(s): %s ...\n"%(len(regions), ",".join(regions)[:3]))
    if data.shape[0]<1:
        logger("[mod_plot][ERROR]  %s row(s) found in %s\n"%(data.shape[0], infn))
        return
    
    # process regions/chromosomes
    ## this can be run in parallel easily, but is it worth the effort really?
    for ref, chrdata, strand in zip(regions, regionsData, strands):
        # define output
        fn = "%s"%ref
        if args.mod: fn = "%s.%s"%(ref, args.mod)
        outfn = os.path.join(outdir, fn)
        # get data frame
        df = bam2calls(bamfiles, chrdata, args.mapq, args.minModProb, MaxPhredProb, can2mods,
                       maxDepth=args.maxDepth, minAlgFrac=args.minAlgFrac, strand=strand)
        #print(df.head()); print(df.columns); print(df.index)
        # get colors for rows and columns
        samples = df.pop("sample")#; print(samples.unique())
        pal = sns.cubehelix_palette(samples.unique().size, light=.9, dark=.1, reverse=True, start=1, rot=-2)
        lut = dict(zip(sorted(samples.unique(), key=lambda x: x.split("_")[-1]), pal))
        row_colors = [samples.map(lut), ]
        read_strand = df.pop("read_strand")#; print(read_strand.unique())
        if len(read_strand.unique())>1:            
            row_colors.append(read_strand.map(read_strand_lut))
        mods = df.columns.get_level_values("mod")
        mods_lut = dict(zip(sorted(mods.unique()), args.colors))
        col_colors = mods.map(mods_lut)
        # https://seaborn.pydata.org/generated/seaborn.clustermap.html
        g = sns.clustermap(df, cmap="Blues", method=args.method, metric=args.metric, figsize=(16, 10), 
                           col_cluster=args.col_cluster, col_colors=col_colors, row_colors=row_colors, 
                           xticklabels=False, yticklabels=False, 
                           cbar_kws={'label': 'Modification probability', 'orientation': 'horizontal'},)
        # add legends
        l1 = g.fig.legend(loc='lower left',bbox_to_anchor=(0.01, 0.8), frameon=True, title="Sample", 
                          handles=[mpatches.Patch(color=c, label=l) for l, c in lut.items()])
        l2 = g.fig.legend(loc='lower right', bbox_to_anchor=(1.0, 0.8), frameon=True, ncol=10, title="Modification", 
                          handles=[mpatches.Patch(color=c, label=l) for l, c in mods_lut.items()])
        if len(read_strand.unique())>1: 
            l3 = g.fig.legend(loc='upper left',bbox_to_anchor=(0.01, 0.8), frameon=True, ncol=1, title="Read strand", 
                              handles=[mpatches.Patch(color=c, label=l) for l, c in read_strand_lut.items()])
        # adjust legend horizontal and top left corner
        g.cax.set_position([.01, .05, .15, .03])#[.8, .96, .2, .03])
        # set title and labels for axes
        g.fig.suptitle("Clustermap for %s %s %s"%(ref, args.method, args.mod)) #g.ax_col_dendrogram.set_title(
        g.ax_heatmap.set_xlabel("Modified positions")
        g.ax_heatmap.set_ylabel("Reads")
        # save
        g.savefig("%s.clustermap.%s.%s"%(outfn, args.method, ext))

        # pca
        pc_op = PCA()
        data_pcs = pc_op.fit_transform(g.data) 
        fig, ax = plt.subplots(1, figsize=(6,5))
        # plot explained variance as a fraction of the total explained variance
        ax.plot(np.arange(1, len(pc_op.explained_variance_)+1),
                pc_op.explained_variance_/pc_op.explained_variance_.sum())
        ax.set_xlabel('Component number')
        ax.set_ylabel('Fraction of explained variance')
        ax.set_title('Scree plot for %s'%ref)
        #fig.tight_layout()
        fig.savefig("%s.scree.%s.%s"%(outfn, args.method, ext)) #show()
        ax.set_xlim((0, 10))
        fig.savefig("%s.scree_10.%s.%s"%(outfn, args.method, ext)) #show()
        # clean-up
        g.fig.clear()
        plt.close()

def collapse_axes(xlab, ylab):
    """Return X/Y labels for unique X positions"""
    _xlab, _ylab, pos = [xlab[0]], [ylab[0]], [0]
    for i in range(1, len(xlab)):
        if xlab[i] != xlab[i-1]:
            _xlab.append(xlab[i])
            _ylab.append(ylab[i])
            pos.append(i)
    return _xlab, _ylab, pos
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    methods = ["single", "ward", "complete", "average", "weighted", "centroid", "median"]
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
    metrics = ["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine",
               "dice", "euclidean", "hamming", "jaccard", "jensenshannon", "kulsinski",
               "mahalanobis", "matching", "minkowski", "rogerstanimoto", "russellrao",
               "seuclidean", "sokalmichener", "sokalsneath", "sqeuclidean", "yule"]
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--input", default="modPhred/mod.gz", help="input file [%(default)s]")
    parser.add_argument("-o", "--outdir", default="modPhred/correlations", help="output dir [%(default)s]")
    parser.add_argument("-b", "--bams", default=[], nargs="*", help="input BAMs [read from mod.gz]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="min mapping quality [%(default)s]")
    parser.add_argument("-d", "--mindepth", default=25, type=int, help="min depth of coverage [%(default)s]")
    parser.add_argument("--maxDepth", default=100, type=int, help="max depth of coverage [%(default)s]")
    parser.add_argument("--minfreq", "--minModFreq", default=0.20, type=float, help="min modification frequency per position [%(default)s]")
    parser.add_argument("--minModProb", default=0.50, type=float, help="min modification probability per base [%(default)s]")
    parser.add_argument("--minAlgFrac", default=0.80, type=float, help="min fraction of read aligned to the region [%(default)s]")
    parser.add_argument("--col_cluster", action="store_true", help="cluster also columns")    
    parser.add_argument("--colors", default="cmybrgwk", help="colors for plotting [%(default)s]")
    #parser.add_argument("-s", "--strand", default=None, choices=["+", "-"], help="select strand [include both]")
    parser.add_argument("--mod", default="", help="filter only 1 modification [analyse all]")
    parser.add_argument("--method", default="ward", choices=methods, help="method used for clustering [%(default)s]")    
    parser.add_argument("--metric", default="euclidean", choices=metrics, help="metric used for clustering [%(default)s]")    
    parser.add_argument("-r", "--regions", nargs="+", default=[], help="regions to process [all chromosomes]")    
    parser.add_argument("-w", "--overwrite", action="store_true", help="overwrite existing output")    
    parser.add_argument("-e", "--ext", default="svg", help="figure format/extension [%(default)s]")
    
    o = parser.parse_args() 
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    if not o.regions:
        logger("Processing entire chromosomes - consider narrowing to certain regions!")
        
    mod_cluster(o.outdir, o.input, o.bams, args=o, ext=o.ext)
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

#!/usr/bin/env python3
desc="""Generate plots based on var.tsv.gz (mod_report.py output). 

More info at: https://github.com/lpryszcz/modPhred

Dependencies: numpy, pandas, matplotlib
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 23/06/2019
"""

import gzip, os, pickle, sys
from datetime import datetime
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
import matplotlib.pyplot as plt
import seaborn as sns
from guppy_encode import VERSION, HEADER, load_info, logger

def get_positions_from_bed(fname):
    """Return list of modified positions"""
    positions = []
    for l in open(fname):
        if l.startswith('#') or not l[:-1]:
            continue
        ldata = l[:-1].split('\t')
        positions.append(":".join((ldata[0], ldata[2])))
    return positions

def plot_venn(outfn, beds, names=[], title=""):
    """Plot venn diagram"""
    import venn
    # select plotting function
    if len(beds)==2: func = venn.venn2
    elif len(beds)==3: func = venn.venn3
    elif len(beds)==4: func = venn.venn4
    elif len(beds)==5: func = venn.venn5
    elif len(beds)==6: func = venn.venn6
    else:
        logger("[ERROR] Please provide between 2 and 6 BED files\n")
        sys.exit(1)
    # use fnames as names if names not given
    if len(names)!=len(beds): names = beds    
    # load positions & plot
    labels = venn.get_labels([get_positions_from_bed(bed) for bed in beds])#, fill=['number', 'logic']
    fig, ax = func(labels, names=names)
    # add title
    if title: plt.title(title)
    # and save or visualise plot
    if outfn: plt.savefig(outfn)
    else: plt.show()
   
def load_bed(fname):
    """Return regions from BED file. If not a file, try to unload region(s) from a string"""
    regions = []
    if os.path.isfile(fname) or os.path.islink(fname):
        for l in open(fname):
            if l.startswith('#') or not l[:-1]:
                continue
            ldata = l[:-1].replace(',','').split('\t')#; print ldata
            if len(ldata) >= 3:
                ref, start, end = ldata[:3]
            else:
                ref, se = ldata[0].split(':')
                start, end = se.split('-')
            start, end = map(int, (start, end))
            regions.append((ref, start, end)) #yield ref, start, end
    else:
        for region in fname.split():
            if not region: continue
            ref, se = region.replace(',','').split(':')
            start, end = se.split('-')
            start, end = map(int, (start, end))
            regions.append((ref, start, end))
    return regions

def plot_scatter(infn, ext="png", logger=logger, data=False, region="", samples=[], 
                 features=["depth", "basecall_accuracy", "mod_frequency", "median_mod_prob"]):
    """Plot scatter using seaborn"""
    # make sure outfn exists
    if not os.path.isfile(infn):
        logger("[mod_plot][ERROR] File %s does not exists! Have you run mod_report.py?\n"%infn)
        sys.exit(1)
    # get outdir
    outdir = os.path.join(os.path.dirname(infn), "plots")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    logger("Saving plots for %s to %s ...\n"%(", ".join(features), outdir))
        
    # parse data
    if isinstance(data, bool):
        logger("Loading %s ...\n"%infn)
        data = pd.read_csv(infn, sep="\t", header=len(HEADER.split('\n'))-2, index_col=False,
                           dtype={"chr": object, "pos": int}) # ADD TO ALL
    # limit by region AND CONSIDER LIMITING COV TO 2-3x median?
    if region:
        logger(" limiting to %s ...\n"%region)
        chrom, s, e = region, 0, 0
        if "-" in region:
            chrom, se = region.split(':')
            s, e = map(int, se.split("-"))
        data = data[data.chr==chrom]
        if e:
            data = data[s<=data.pos<=e]
    if data.shape[0]<1:
        logger("[mod_plot][ERROR]  %s row(s) found in %s\n"%(data.shape[0], infn))
        return
    # rename .bam columns to basename and split it at _ with \n
    data.columns = [os.path.basename(c).replace("_","\n") if ".bam" in c else c for c in data.columns]#; data.head()
    # plot features
    for feature in features:
        # match feature by replacing _ with \n
        cols = list(filter(lambda c: feature.replace("_","\n") in c, data.columns))#; print(cols)
        # limit by sample
        if samples:
            scols = set(c for c in cols for s in samples if s.replace("_","\n") in c)
            cols = list(sorted(scols))#; print(cols)
        # plot palette="husl", 
        g = sns.pairplot(data, vars=cols, height=4, hue='mod', diag_kind='kde',
                         plot_kws={'alpha': 0.1, 's': 3, }) #'edgecolor': 'k'
        outfn = os.path.join(outdir, "%s.%s"%(feature, ext))
        # add figure title
        g.fig.suptitle("%s\n%s"%(feature, infn), size=16)
        g.fig.subplots_adjust(top=.90, right=0.95)
        # make legend in top right corner and increase marker size
        g._legend.set_bbox_to_anchor((0.10, 0.95))
        g._legend.set_title("") #"mods:"
        for lh in g._legend.legendHandles: lh._sizes = [20] #lh.set_alpha(1)
        # set axes limits 0-1
        if feature!="depth":
            for r in g.axes: 
                for ax in r:
                    ax.set_xlim((0, 1))
                    ax.set_ylim((0, 1))
        # save
        g.fig.savefig(outfn)
        
def plot_regions(infn, bed, ext="svg", logger=logger, data=False, colors="brcmyg"): 
    """Generate frequency plots for given regions

    If data is given, it won't be loaded again.
    """
    # make sure outfn exists
    if not os.path.isfile(infn):
        logger("[mod_plot][ERROR] File %s does not exists! Have you run mod_report.py?\n"%infn)
        sys.exit(1)
    # get outdir
    outdir = os.path.join(os.path.dirname(infn), "plots")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # load regions to plot
    regions = load_bed(bed)#; print(regions)
    logger("Saving plots for %s region(s) to %s ...\n"%(len(regions), outdir))
        
    # parse data
    if isinstance(data, bool):
        logger("Loading %s ...\n"%infn)
        data = pd.read_csv(infn, sep="\t", header=len(HEADER.split('\n'))-2, index_col=False,
                           dtype={"chr": object, "pos": int})
    if data.shape[0]<1:
        logger("[mod_plot][ERROR]  %s row(s) found in %s\n"%(data.shape[0], infn))
        return
    # get uniq mods
    mods = data["mod"].unique()
    metrics = ['depth', 'basecall_accuracy', 'mod_frequency', 'median_mod_prob']
    sample_names = [get_sample_name(n) for n in data.columns if n.endswith(metrics[0])]
    metric2cols = {m: [c for c in data.columns if c.endswith(m)] for m in ['mod_frequency', ]}#; print(metric2cols)
    logger(" %s samples and %s modifications: %s\n"%(len(sample_names), len(mods), ", ".join(mods)))
    metric = 'mod_frequency'
    # plot regions
    for ref, s, e in regions:
        #df = data[(data.chr==ref)&(data.pos>=s)&(data.pos<=e)]
        df = data[np.all((data.chr==ref, data.pos>=s, data.pos<=e), axis=0)]
        if df.shape[0]<1:
            logger("[mod_plot][ERROR] No modifications in %s:%s-%s\n"%(ref, s, e))
            continue
        mods = df["mod"].unique()
        logger(" %s:%s-%s with %s modifications: %s\n"%(ref, s, e, len(mods), ", ".join(mods)))
        #return
        fig, axes = plt.subplots(nrows=len(sample_names), ncols=1, sharex="all", sharey="all", 
                                 figsize=(20, 2+1*len(sample_names))) #20, 12 for 2 samples
        fig.suptitle("%s:%s-%s"%(ref, s, e), fontsize=12)
        labels = []
        for strand, norm in zip("+-", (1, -1)):
            for ax, col, name in zip(axes, metric2cols[metric], sample_names):
                for color, mod in zip(colors, mods):
                    selection = (df.strand==strand)&(df["mod"]==mod)
                    ax.bar(df[selection].pos, norm*df[selection][col], color=color, label=mod)
                    if strand:
                        ax.set_title(name.replace('\n', '_')) #col)
                        ax.set_ylabel("%s\n[on +/- strand]"%metric)

        # set limits
        ax.set_xlim(s, e+1)
        ax.set_ylim(-1, 1)
        ax.set_xlabel("%s position"%ref)
        #https://stackoverflow.com/a/13589144/632242
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        #plt.legend(by_label.values(), by_label.keys())
        fig.legend(handles=by_label.values(), labels=by_label.keys())
        #fig.tight_layout()        
        #plt.show()
        outfn = os.path.join(outdir, "%s:%s-%s.%s"%(ref, s, e, ext))
        fig.savefig(outfn)

def plot_flipflop(sig, trace, read, seq, s=0, n=20, BASE="ACGTZ", colours={'A':'green', 'C':'blue', 'G':'orange', 'T':'red', 'Z':'purple', 'N': 'grey'}):
    """Plot flipflop trace from Fast5"""
    pos = np.argwhere(move==1).flatten()
    down_sample_factor = round(len(sig) / float(len(trace)))
    # get subset of trace & normalise to 0-1
    trace = trace[pos[s]:pos[s+n]] / 255.
    # get number of bases
    nbase = trace.shape[1] // 2
    # start figure
    fig = plt.figure(figsize=(20, 5))
    # plot squiggle
    ax1 = fig.add_subplot(211)
    plt.title("%s (bases: %s-%s)"%(read, s, s+n))
    plt.ylabel('Normalised signal')
    plt.plot(np.arange(pos[s]*down_sample_factor, pos[s+n]*down_sample_factor), sig[pos[s]*down_sample_factor:pos[s+n]*down_sample_factor], color='grey')
    # plot flipflop
    fig.add_subplot(212, sharex=ax1)
    plt.xlabel('time (samples)')
    plt.ylabel('State probability')
    x2 = down_sample_factor * np.arange(pos[s], pos[s+n])
    for i in range(nbase):
        plt.fill_between(x2, trace[:, i], color=colours[BASE[i]], alpha=0.3)
        plt.fill_between(x2, trace[:, i + nbase], color=colours[BASE[i]], alpha=0.3)
        plt.plot(x2, trace[:, i], color=colours[BASE[i]], label=BASE[i])
        plt.plot(x2, trace[:, i + nbase], color=colours[BASE[i]], linestyle='dashed') 
    # add bases
    for pi, p in enumerate(pos[s:s+n], s):
        plt.text(p*down_sample_factor, -0.1, seq[pi])
    # add legend
    plt.legend()
    plt.grid()
    plt.show()

def violin_plot(data, title="", axis="", names=[]):
    """Generate violin plot"""
    if axis: 
        ax=axis
    else: 
        ax = plt.subplot() 
    # violin plot
    ax.violinplot(data, range(len(data)), points=20, widths=0.7, bw_method=0.5, showmeans=True, showextrema=True, showmedians=True)
    ax.set_xticks(range(len(data)))
    if axis: 
        ax.set_xticklabels([" " for x in range(len(data))])
    elif names:
        ax.set_xticklabels(names)        
    if title:
        ax.set_title(title)
    if not axis: plt.show()
    return ax

def mod_plot(infn, ext="svg", logger=logger, data=False, colors="brcmyg"): #png
    """Generate violin plots

    If data is given, it won't be loaded again.
    """
    # make sure outfn exists
    if not os.path.isfile(infn):
        logger("[mod_plot][ERROR] File %s does not exists! Have you run mod_report.py?\n"%infn)
        sys.exit(1)
        
    # parse data
    if isinstance(data, bool):
        logger("Loading %s ...\n"%infn)
        data = pd.read_csv(infn, sep="\t", header=len(HEADER.split('\n'))-2, index_col=False)
    if data.shape[0]<10:
        logger("[mod_plot][ERROR]  %s row(s) found in %s\n"%(data.shape[0], infn))
        return
    # plot
    bases = data["mod"].unique()#; print(bases)
    metrics = ['depth', 'basecall_accuracy', 'mod_frequency', 'median_mod_prob']
    metrics_names = ["Number of reads", "Agreement with reference",
                     "Frequency of modification", "Median modification probability"]
    sample_names = [get_sample_name(n) for n in data.columns if n.endswith(metrics[0])]
    fig, axes = plt.subplots(nrows=len(metrics), ncols=len(bases), sharex="col", sharey="row", 
                             figsize=(1.5*len(bases)*len(sample_names), 2+3*len(metrics))) #6, 20
    fig.suptitle(infn, fontsize=12)
    nans = [float('nan'), float('nan')]
    # get max median depth
    maxYdepth = 0
    for bi, b in enumerate(bases):
        # get mask for only median_mod_prob
        cols = list(filter(lambda x: x.endswith(metrics[-1]), data.columns))
        _data = data[data["mod"]==b].loc[:, cols].to_numpy()
        # mask nan before plotting https://stackoverflow.com/a/44306965/632242
        mask = ~np.isnan(_data)
        for mi, m in enumerate(metrics):
            cols = list(filter(lambda x: x.endswith(m), data.columns))
            ax = axes[mi, bi] if len(bases)>1 else axes[mi]
            _data = data[data["mod"]==b].loc[:, cols].to_numpy()#; print(bi, b, mi, m, _data.shape)
            #if _data.sum():
            a = ax.violinplot([d[m] if d[m].any() else nans for d, m in zip(_data.T, mask.T)], points=20, widths=0.7,
                              bw_method=0.5, showextrema=True, showmedians=True) #showmeans=True,
            # color samples differently
            for pci, pc in enumerate(a['bodies']): pc.set_facecolor(colors[pci%len(colors)])
            #pc.set_edgecolor('black')
            #pc.set_alpha(1)            
            ax.set_xticks(range(1, len(cols)+1))
            ax.set_xticklabels([" " for x in range(len(cols))])
            if not mi: 
                ax.set_title("%s\n%s positions"%(b, data[data["mod"]==b].shape[0]))
                # set depth Y range as 2*median of depth for A
                if 2*np.nanmedian(_data, axis=0).max()>maxYdepth:
                    maxYdepth = 2*np.nanmedian(_data, axis=0).max()#; print(np.nanmean(_data, axis=0)); print(a['cmedians'])
                    ax.set_ylim(0, maxYdepth)
            elif mi in (1, 3): ax.set_ylim(0.5, 1)
            #elif mi==2:  ax.set_ylim(0, 0.5)
            else: ax.set_ylim(0, 1)
            if not bi: 
                ax.set_ylabel(metrics_names[mi])
            if mi+1 == len(metrics):
                ax.set_xticklabels(sample_names)
            ax.grid(axis="y", which="both")
    #fig.show()
    fig.savefig(infn+".%s"%ext)

def mod_plot_bases(infn, ext="svg", logger=logger, data=False): 
    """Generate violin plots

    If data is given, it won't be loaded again.
    """
    # make sure outfn exists
    if not os.path.isfile(infn):
        logger("[mod_plot][ERROR] File %s does not exists! Have you run mod_report.py?\n"%infn)
        sys.exit(1)
        
    # parse data
    if isinstance(data, bool):
        logger("Loading %s ...\n"%infn)
        data = pd.read_csv(infn, sep="\t", header=len(HEADER.split('\n'))-2, index_col=False)
    if data.shape[0]<10:
        logger("[mod_plot][ERROR]  %s row(s) found in %s\n"%(data.shape[0], infn))
        return
    # plot
    bases = 'ACGT'
    metrics = ['depth', 'basecall_accuracy', 'mod_frequency', 'median_mod_prob']
    sample_names = [get_sample_name(n) for n in data.columns if n.endswith(metrics[0])]
    fig, axes = plt.subplots(nrows=len(metrics), ncols=len(bases), sharex="col", sharey="row", 
                             figsize=(2+1.5*len(bases)*len(sample_names), 5*len(metrics))) #6, 20
    fig.suptitle(infn, fontsize=12)
    nans = [float('nan'), float('nan')]
    # get max median depth
    maxYdepth = 0
    for bi, b in enumerate(bases):
        # get mask for only median_mod_prob
        cols = list(filter(lambda x: x.endswith(metrics[-1]), data.columns))
        _data = data[data.ref_base==b].loc[:, cols].to_numpy()
        # mask nan before plotting https://stackoverflow.com/a/44306965/632242
        mask = ~np.isnan(_data)
        for mi, m in enumerate(metrics):
            cols = list(filter(lambda x: x.endswith(m), data.columns))
            ax = axes[mi, bi]
            _data = data[data.ref_base==b].loc[:, cols].to_numpy()#; print(bi, b, mi, m, _data.shape)
            #if _data.sum():
            a = ax.violinplot([d[m] if d[m].any() else nans for d, m in zip(_data.T, mask.T)], points=20, widths=0.7,
                              bw_method=0.5, showmeans=True, showextrema=True, showmedians=True)
            ax.set_xticks(range(1, len(cols)+1))
            ax.set_xticklabels([" " for x in range(len(cols))])
            if not mi: 
                ax.set_title("%s (%s positions)"%(b, data[data.ref_base==b].shape[0]))
                # set depth Y range as 2*median of depth for A
                if 2*np.nanmedian(_data, axis=0).max()>maxYdepth:
                    maxYdepth = 2*np.nanmedian(_data, axis=0).max()#; print(np.nanmean(_data, axis=0)); print(a['cmedians'])
                    ax.set_ylim(0, maxYdepth)
            else:
                ax.set_ylim(0, 1)
            if not bi: 
                ax.set_ylabel(m)
            if mi+1 == len(metrics):
                ax.set_xticklabels(sample_names)
    ax.set_ylim(0.5, 1)
    #fig.show()
    fig.savefig(infn+".%s"%ext)

def get_sample_name(col_name):
    """Return sample name from column name"""
    # modPhred/m6A/minimap2/RNA081120181_unmodified.bam depth > RNA081120181_unmodified.bam
    # modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145.bam > MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145.bam
    fname = os.path.basename(col_name.split()[0]) #
    # RNA081120181_unmodified.bam > unmod
    # RNAAA023484_wt1.bam > wt1
    # RNA081120181_unmodified.bam > unmod
    if fname.startswith(("RNA", "cDNA", "DNA")):
        return fname[:-4].split("_")[-1][:5]
    else:
        return "\n".join(fname[:-4].split("_"))#[-3:]
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-e", "--ext", default="svg", help="figure format/extension [%(default)s]")
    qcargs = parser.add_argument_group('Options for plotting QC plots')
    qcargs.add_argument("-i", "--input", default="modPhred/mod.gz", help="input file [%(default)s]")
    #parser.add_argument("-d", "--depthFrac", default=0.1, type=float,
    #                    help="samples with depth of coverage lower than depthFrac* median will be skipped [%(default)s]")
    region = parser.add_argument_group('Options for plotting frequency in regions')
    region.add_argument("-b", "--bed", help="BED file with regions to analyse [optionally]")
    scatter = parser.add_argument_group('Options for scatter plots (pairplots)')
    scatter.add_argument("--scatter", action='store_true', help="plot scatter plot of multiple features")
    scatter.add_argument("-r", "--region", default="", help="include only mods from given region [all]")
    scatter.add_argument("-s", "--samples", default=[], nargs="+", help="include only samples containing str [all]")
    vennargs = parser.add_argument_group('Options for plotting Venn diagrams')
    vennargs.add_argument("--venn", nargs="+", default=[], help="BED files to plot Venn diagram [optionally, up to 6 files]")
    vennargs.add_argument("-n", "--names", nargs="*", default=[], help="names to plot on Venn [will use file names]")
    vennargs.add_argument("-o", "--out", default=None, help="save Venn as file [render]")
    vennargs.add_argument("-t", "--title", default=None, help="figure title")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    if o.bed:
        plot_regions(o.input, o.bed, ext=o.ext)
    elif o.scatter:
        plot_scatter(o.input, ext=o.ext, region=o.region, samples=o.samples)
    elif o.venn:
        plot_venn(o.out, o.venn, o.names, o.title)
    else:
        # get plots
        mod_plot(o.input, ext=o.ext)#, o.depthFrac)

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

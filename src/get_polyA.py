#!/usr/bin/env python3
desc="""Basecall Fast5 (or read baescalled Fast5), align & and estimate polyA tail length.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 26/08/2021
"""

import csv, gzip, inspect, mappy, numpy as np, os, pysam, sys, tempfile
from datetime import datetime
from pathlib import Path
#from scipy.signal import savgol_filter, find_peaks
from multiprocessing import Pool
from basecall import init_args, start_guppy_server, get_basecall_client, basecall_and_align, get_sam_header
#from resquiggle import resquiggle_reads
from common import VERSION, logger
import numba

old_settings = np.seterr(divide='ignore') #under='ignore')

def get_consecutive(data, stepsize=1):
    """Return consecutive windows allowing given max. step size"""
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)

@numba.jit
def compute_running_mean_diffs(sig, n=7, w=50):
    """Compute average difference between n neighboring window means each of size w"""
    moving_average = np.cumsum(sig) # 19.9 µs
    moving_average[w:] = (moving_average[w:]-moving_average[:-w])
    moving_average = moving_average[w-1:] / w
    # extract moving window averages at n offsets
    offsets = [moving_average[int(w*offset):int(-w*(n-offset-1))]
               for offset in range(n-1)] + [moving_average[int(w*(n-1)):], ]
    # compute difference between all pairwise offset
    diffs = [np.abs(offsets[i]-offsets[j]) for i in range(n) for j in range(i+1, n)]
    # compute average over offset differences at each valid position
    diff_sums = diffs[0] #.copy()
    for diff_i in diffs[1:]: diff_sums += diff_i
    return diff_sums / len(diffs)

def plot(sig, sig_step, s, e, title="", first_base=None):
    """Plot read signal (sig) with changes in signal (sig_step) 
    and mark polyA tail withing s:e boundaries. 
    Plot signal up to first_base position. 
    """
    fig, ax = plt.subplots(figsize=(15, 3))
    ax.plot(sig[:30000])
    ax.axvspan(s, e, color="red", alpha=.3, label="polyA")
    if first_base:
        ax.axvspan(first_base, len(sig), color="cyan", alpha=.3, label="aliged_read")
    ax1 = ax.twinx()
    ax1.plot(sig_step, c="orange")
    ax.legend()
    pA = sig[s:e]
    if title: ax.set_title(title)

def get_polyA_boundary(sig_all, md, read_pos, read_start_trim, rname, plot_read=False, 
                       n=3, window=50, threshold=0.4, include_bases=20):
    """Return polyA start and end in signal space.
    Here, signal should be already trimmed to the first aligned base
    """
    s = e = 0
    flt = "FAIL"
    sig = sig_all[:read_pos[read_start_trim+include_bases]] # 219 ns  
    sig_pA = (sig+md['daq_offset'])*md['daq_scaling'] # 8.32 µs; 28.6 µs for the whole read
    sig_norm = (sig_pA-md["median"])/md["med_abs_dev"] # 7.41 µs; 226.3 µs for the whole read
    # compute running mean differences
    rmd = compute_running_mean_diffs(sig_norm, n, window) # 50µs with numba; 366 µs vs 900µs for fftconvolve
    stalls = []
    # get start and end of stalls
    for el in get_consecutive(np.argwhere(rmd<threshold)[:, 0], window):
        if not len(el): continue
        _s, _e = el[[0, -1]]
        # with at least n*window size
        if _e-_s>n*window:
            # and extend end by n*window
            stalls.append([_s, _e+n*window])
    if len(stalls): 
        merged_stalls = [stalls[0][:], ]
        for i in range(1, len(stalls)):
            # join overlapping stalls - here you may add buffer of window size
            if stalls[i][0]+window<merged_stalls[-1][1]: merged_stalls[-1][1] = stalls[i][1]
            else: merged_stalls.append(stalls[i])
        # filter stalls by expected normalised signal
        stalls = [(s, e) for s, e in merged_stalls if 0<sig_norm[s:e].mean()<3
                  and s>read_pos[read_start_trim]/10 # polyA start not before first_base/10
                  and 3*(e-s)>abs(read_pos[read_start_trim]-e) # polyA distance to first base < 3*polyA length
                  and e!=read_pos[read_start_trim+include_bases] # polyA ending at include_bases
                  ]
    if len(stalls):
        s, e = sorted(stalls, key=lambda x: x[1]-x[0])[-1]
        flt = "PASS" if 0.1*sig[s:e].mean()>sig[s:e].std() else "FAIL"
    # catch very short tails
    if e-s<3: flt="FAIL"
    # plot if requested
    if plot_read:
        print(rname, s, e, flt, np.mean(sig_norm[s:e]))
        plot(sig_norm, rmd, s, e, "%s: %s"%(rname, flt), read_pos[read_start_trim])
    return s, e, flt

methods = ("polyA", "polyA_flt", "mean_flt", "median")
def get_polyA_len(a, sig, move, md, rna=True, perc=5,
                  op_consuming_reference=set((0, 2, 7, 8))):
    """Return filter (FAIL or PASS) and tuple of estimation of polyA lenght using:
    - mean
    - median
    - median of only A bases
    - find_peaks

    It returns 0 for all methods if the estimation fails. 
    """
    # remove outliers from signal - is it really needed here?
    #sig[(sig<0)|(sig>2000)] = np.mean(sig)
    # get positions of read bases in signal space
    read_pos = md['trimmed_samples']+np.argwhere(move==1)[:, 0]*md['model_stride']
    # get number of bases trimmed from end and start of the read
    ## this is for RNA - meaning read is reversed (3'>5')
    ## be careful if you want to analyse cDNA!
    read_start_trim = a.cigar[-1][1] if a.cigar[-1][0]==4 else 0
    read_end_trim = a.cigar[0][1] if a.cigar[0][0]==4 else 0
    if a.is_reverse: 
        read_start_trim, read_end_trim = read_end_trim, read_start_trim
    # get polyA tail boundiaries
    s, e, flt = get_polyA_boundary(sig, md, read_pos, read_start_trim, a.reference_name)
    # estimated aligned reference lenght - this should be intro free
    rlen = sum(c for op, c in a.cigar if op in op_consuming_reference)
    # get polyA length using simple ratio of pA signal over total signal x number of bases
    sig_len = read_pos[len(a.seq)-read_end_trim-1] - read_pos[read_start_trim]
    polyA_len = round((e-s)/sig_len*rlen, 2)
    # using dt mean flt
    dt = read_pos[1:]-read_pos[:-1]; dt
    dt_clean = dt[(dt>=np.percentile(dt, perc))&(dt<=np.percentile(dt, 100-perc))]    
    mean_flt = round((e-s)/dt_clean.mean(), 2)
    # and median
    median = round((e-s)/np.median(dt_clean), 2)
    # get rid of likely stalled estimates
    dt_flt = dt[(dt<=np.percentile(dt, perc))|(dt>=np.percentile(dt, 100-perc))]
    polyA_len_flt = round((e-s)/(sig_len-dt_flt.sum())*(rlen-len(dt_flt)), 2)
    return flt, (polyA_len, polyA_len_flt, mean_flt, median)

def get_polyA_from_fast5(args):
    """Process individual Fast5 files"""
    ref, fn, conf, rna = args
    header_dict = get_sam_header(ref)#, coord_sorted=True)
    header = pysam.AlignmentHeader.from_dict(header_dict)
    i = 0
    flt = "PASS"
    data = []
    for a, sig, move, modbaseprob, md in basecall_and_align(fn, header, rna, conf):
        if a.is_unmapped: continue
        i += 1
        #if i>200: break
        if not i%100: sys.stderr.write(" %s  \r"%i)
        # allow secondary if primary failed
        if flt=="PASS" and a.is_supplementary or a.is_secondary: continue
        # add polyA estimation
        flt, polyA_est = get_polyA_len(a, sig, move, md, rna)
        data.append((a.qname, a.reference_name, flt, *polyA_est))
    return data

def get_polyA(indirs, fasta, threads=1, rna=True, config="dna_r9.4.1_450bps_pcr_fast",
              host="localhost", port=5555, recursive=False, device="cuda:0"):
    """Process multiple directories from Fast5 files"""
    # start guppy server if needed
    conf = (config, host, port)
    if host:
        guppy_proc, host, port = start_guppy_server(host, config, port, device)
        # try connecting to guppy server first
        client, _get_read_data = get_basecall_client(*conf)
    else:
        guppy_proc = None
        logger("We'll use basecall information from Fast5 files...")
    # load reference for mappy
    logger("Loading reference index from %s..."%fasta)
    aligner = mappy.Aligner(fasta, preset="spliced" if rna else "map-ont", k=13)
    # start pool of workers
    ## we can't use Queue here, since guppy results (read2data) will be ofter too large!
    # it's important to initialise the pool with aligner object as it can't be pickled
    p = Pool(threads, initializer=init_args, initargs=(aligner, )) # maxtasksperchild=1
    # process directories
    logger("Processing Fast5 files in %s directories...\n"%len(indirs))
    for indir in indirs:
        outfn = "%s.modPhred.polyA.tsv.gz"%indir
        if os.path.isfile(outfn):
            logger(" %s exists"%outfn)            
            continue
        # write header
        out = gzip.open(outfn, "wt")
        tsv = csv.writer(out, delimiter='\t')
        header = ["read_id", "ref", "filter", *methods]
        tsv.writerow(header)
        # get Fast5 files
        if recursive:
            fnames = sorted(map(str, Path(indir).rglob('*.fast5')))
        else:
            fnames = sorted(map(str, Path(indir).glob('*.fast5')))
        logger(" %s with %s Fast5 file(s)...\n"%(indir, len(fnames)))    
        # basecall individual fast5 on-the-fly
        args = [(fasta, fn, conf, rna) for fn in fnames]# if not os.path.isfile("%s.bam"%fn)]
        for ii, data in enumerate(p.imap(get_polyA_from_fast5, args), 1):
            sys.stderr.write(" %s / %s \r"%(ii, len(args)))
            if data:
                tsv.writerows(data)
        out.close()
    # close pool
    p.close()
    # and guppy basecall server if it was started by this process
    if guppy_proc: guppy_proc.terminate()
    logger("Done!")
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--indirs", nargs="+", help="input directory with Fast5 files")
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-t", "--threads", default=6, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-c", "--config", default="rna_r9.4.1_70bps_ivt_fast.cfg", help="guppy model [%(default)s]")
    parser.add_argument("--host", "--guppy_basecall_server", default="",
                        help="guppy server hostname or path to guppy_basecall_server binary [use basecall information from Fast5]")
    parser.add_argument("-p", "--port", default=5555, type=int,
                        help="guppy server port (this is ignored if binary is provided) [%(default)s]")
    parser.add_argument("-d", "--device", default="cuda:0", help="CUDA device to use [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    get_polyA(o.indirs, o.fasta, o.threads, o.rna,
              o.config, o.host, o.port, o.recursive, o.device)
        
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


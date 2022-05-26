#!/usr/bin/env python3
desc="""Convert basecalled Fast5 with modifications annotated by guppy v3.1.5+
to FastQ with modification probabilities encoded as FastQ qualities.

More info at: https://github.com/lpryszcz/modPhred

Dependencies: h5py

TO DO:
- drop low quality reads
- ignore low quality bases
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 20/06/2019
"""

import glob, gzip, h5py, os, pickle, resource, sys, subprocess, time
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path
import numpy as np

VERSION = '1.0c'
# set max 3 modifications per each base
MaxModsPerBase = 3
MaxProb = 256
QUALS = "".join(map(chr, range(33, 127)))
# it's only DNA as in SAM U should be A
base2complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

HEADER = """###
# Welcome to modPhred (ver. %s)!
# 
# Executed with: %s
#
# For each bam file 4 values are stored for every position:
# - depth of coverage (only positions with >=%s X in at least one sample are reported)
# - accuracy of basecalling (fraction of reads having same base as reference, ignoring indels)
# - frequency of modification (fraction of reads with modification above given threshold)
# - median modification probability of modified bases (0-1 scaled). 
#
# If you have any questions, suggestions or want to report bugs,
# please use https://github.com/novoalab/modPhred/issues.
# 
# Let's begin the fun-time with Nanopore modifications...
###
chr\tpos\tref_base\tstrand\tmod\t%s
"""

# check if all executables exists & in correct versions
def byte2str(e):
    """Return str representation of byte object"""
    return e.decode("utf-8")

def _check_executable(cmd):
    """Check if executable exists."""
    p = subprocess.Popen("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return "".join(map(byte2str, p.stdout.readlines()))

def _check_dependencies(dependencies):
    """Return error if wrong software version"""
    warning = 0
    # check dependencies
    info = "[WARNING] Old version of %s: %s. Update to version %s+!\n"
    for cmd, version in dependencies.items():
        out = _check_executable(cmd)
        if "not found" in out: 
            warning = 1
            sys.stderr.write("[ERROR] %s\n"%out)
        elif version:
            p = subprocess.Popen([cmd, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out = byte2str(p.stdout.readlines()[0])
            curver = out.split()[-1].split('-')[0] # 2.16-r922 for minimap2
            try:
                curver = float(curver)
                if curver<version:
                    warning = 1
                    sys.stderr.write(info%(cmd, curver, version))
            except:
                warning = 1
                sys.stderr.write("[WARNING] Problem checking %s version: %s\n"%(cmd, out))
                
    message = "Make sure you have installed all dependencies from https://github.com/lpryszcz/modPhred#prerequisites !"
    if warning:
        sys.stderr.write("\n%s\n\n"%message)
        sys.exit(1)
    
dependencies = {'minimap2': 2.16, 'samtools': 1.3}
#_check_dependencies(dependencies)

def get_MaxPhredProb(MaxModsPerBase, QUALS=QUALS):
    MaxPhredProb = int(len(QUALS)/MaxModsPerBase)
    return MaxPhredProb

def dump_info(indir, alphabet, symbol2modbase, canonical2mods, base2positions, fast5, fast5mod, 
              MaxModsPerBase, MaxPhredProb, fn="modPhred.pkl"):
    """Dump Fast5 info to pickle."""
    data = {"version": VERSION,
            "MaxModsPerBase": MaxModsPerBase,
            "MaxProb": MaxProb,
            "QUALS": QUALS,
            "MaxPhredProb": MaxPhredProb,
            "alphabet": alphabet,
            "symbol2modbase": symbol2modbase,
            "canonical2mods": canonical2mods,
            "base2positions": base2positions,
            "fast5": fast5,
            "fast5mod": fast5mod, 
            }
    # dump using pickle
    fname = os.path.join(indir, fn)
    with open(fname, "wb") as out:
        pickle.dump(data, out, protocol=pickle.HIGHEST_PROTOCOL)

def load_info(indir, recursive=False, fn="modPhred.pkl"):
    """Load and return Fast5 info from pickle. Return empty dict if no pickle"""
    data = {} ## ADD RECURSIVE INDIR HANDLING HERE!!
    fname = os.path.join(indir, fn)
    if os.path.isfile(fname):
        data = pickle.load(open(fname, "rb"))
    return data

def get_moltype(output_alphabet):
    """Return moltype and underlying bases for given alphabet"""
    if "T" in output_alphabet:
        bases = "ACGT"
        moltype = "DNA"
    elif "U" in output_alphabet:
        bases = "ACGU"
        moltype = "RNA"
    else:
        logger("[ERROR] Cannot guess if Fast5 contains DNA or RNA (%s)!"%"".join(output_alphabet))
        sys.exit(1)
    return moltype, bases
    
def is_rna(indir):
    """Return True if RNA sample"""
    data = load_info(indir)
    if "U" in data["alphabet"]:
        return True
    return False
    
def get_alphabet(output_alphabet, mods, canonical_bases="ACGTU", force_rna=0):
    """Return correctly ordered alphabet of bases with modified bases
    and 2 dictionaries: symbol2modbase & canonical2mods.
    """
    mods = mods.split() # list of base_mod_long_names ie 6mA 5mC
    output_alphabet = list(output_alphabet) # AYCZGT
    # get ordered alphabet and dictionaries
    symbol2modbase, canonical2mods = {}, {}
    alphabet, _mods = [], []
    # decide if DNA or RNA # ABCGH4TW
    if force_rna: output_alphabet[output_alphabet.index("T")] = "U"
    bases = "".join(b for b in output_alphabet if b in canonical_bases)
    canonical2mods = {b: [] for b in bases}
    idx = -1 
    for b in output_alphabet:
        alphabet.append(b)
        if b in bases:
            idx += 1
        else:
            if b.isdigit():
                canonical2mods[bases[idx+1]].append(b)
            else:
                canonical2mods[bases[idx]].append(b)
            symbol2modbase[b] = mods.pop(0)
    # get base2positions
    base2positions = {}
    idx = 0
    for b in bases:
        base2positions[b] = [idx, ]
        idx += 1
        for mb in canonical2mods[b]:
            base2positions[b].append(idx)
            idx += 1
    #print(bases, alphabet, symbol2modbase, canonical2mods, base2positions)
    return alphabet, symbol2modbase, canonical2mods, base2positions

def get_phredmodprobs(seq, modbaseprobNorm, mods2count, base2positions, canonical2mods, MaxPhredProb):
    """Return PHRED scaled probability of given base being modified"""
    # prepare phred str as numpy array
    phredmodprobs = np.empty(len(seq), dtype='|S1') #array(([QUALS[prob+MaxPhredProb*idx]]*len(seq))
    phredmodprobs[:] = QUALS[0]
    # get positions of each base that can be modified
    base2idx = {b: [] for b in base2positions if len(base2positions[b])>1}
    for ii, b in enumerate(seq): # this can be likely faster with mask or just seq=="A"
        if b in base2idx:
            base2idx[b].append(ii)
    # calculate PHRED scaled modification probability for each base that can be modified
    for b, ii in base2idx.items():
        # only calculate indices max if more than 1 modification per base
        if len(base2positions[b])>2:
            s, e = base2positions[b][1], base2positions[b][-1]+1
            #print(b, modbaseprobNorm[ii].shape, s, e, len(ii))#; print(modbaseprobNorm[ii,s-1:e])
            probs = modbaseprobNorm[ii, s:e]
            indices = np.argmax(probs, axis=1)
            probs = np.max(probs, axis=1)
        else:
            probs = modbaseprobNorm[ii, base2positions[b][1]] 
            indices = np.zeros(probs.shape[0], dtype='int')
        #print(probs, indices)
        # update PHRED
        phredmodprobs[ii] = [QUALS[v] for v in probs+indices*MaxPhredProb]
        # calculate stats
        for idx in range(len(base2positions[b])-1):
            mods2count[canonical2mods[b][idx]] += np.sum(probs[indices==idx]>=0.5*MaxPhredProb)
    return phredmodprobs.tobytes().decode(), mods2count

def get_latest_basecalling(h5, rname):
    """Return latest basecalling from Fast5 file"""
    basecalls = list(filter(lambda x: x.startswith("Basecall_1D_"), h5[rname]['Analyses'].keys()))
    #print(h5[rname]['Analyses'].keys(), basecalls)
    return basecalls[-1]

def worker(args):
    """mod_encode worker using data stored by guppy3"""
    fn, ofn, basecall_group, MaxPhredProb = args
    basecount = warns = 0
    data, rname = [], ""
    alphabet, symbol2modbase, canonical2mods, base2positions, mods2count = '', {}, {}, {}, {}
    # open out file with gzip compression
    outfn = ofn+".fastm.gz"
    outfq = ofn+".fastq.gz_"
    # The default mode is "rb", and the default compresslevel is 9.
    out = gzip.open(outfn, "wt")
    out2 = gzip.open(outfq, "wt")
    # process entries in fast5
    h5 = h5py.File(fn, 'r')
    for i, rname in enumerate(h5):
        #if i>100: break
        # try to load base modification probabilities
        try:                
            if not basecall_group:
                basecall_group = get_latest_basecalling(h5, rname)
            modbaseprob = h5["%s/Analyses/%s/BaseCalled_template/ModBaseProbs"%(rname, basecall_group)]
        except:
            warning("[mod_encode][WARNING] Can't get base modification probabilities for %s:%s!"%(fn, rname))
            warns += 1
            break #continue
        # prepare data storage if not already prepared
        if not data:
            data = [[] for i in range(modbaseprob.shape[1])]
            mods = modbaseprob.attrs['modified_base_long_names'].decode()
            output_alphabet = modbaseprob.attrs['output_alphabet'].tobytes().decode()
            alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(output_alphabet, mods)
            rna = True if "U" in alphabet else False 
            mods2count = {m: 0 for mods in canonical2mods.values() for m in mods}
        # get mod probabilities and normalise to MaxPhredProb
        modbaseprob = np.array(modbaseprob, dtype='float')
        modbaseprobNorm = np.array(modbaseprob * MaxPhredProb / MaxProb, dtype='uint8')
        # reverse only if RNA
        if rna: modbaseprobNorm = modbaseprobNorm[::-1]
        # get read name and sequence as list
        fastq = h5["%s/Analyses/%s/BaseCalled_template/Fastq"%(rname, basecall_group)]
        # unload fastq entry from hdf5 as string, decode and split to list of individual lines
        fastq_elements = fastq[()].tobytes().decode().split('\n')
        read_name = fastq_elements[0].split()[0][1:]
        seq = fastq_elements[1]
        qual = fastq_elements[3]
        basecount += len(seq)
        # get modprobs as qualities
        phredmodprobs, mod2count = get_phredmodprobs(seq, modbaseprobNorm, mods2count, base2positions, canonical2mods, MaxPhredProb)
        # store FastQ and FastM
        out2.write("@%s\n%s\n+\n%s\n"%(read_name, seq, qual))
        out.write("@%s\n%s\n+\n%s\n"%(read_name, seq, phredmodprobs))
    # report number of bases
    if rname:
        name = "/".join(fn.split("/")[-4:])
    # if warnings skip file entirely
    if warns:
        return fn, 0, mods2count, alphabet, symbol2modbase, canonical2mods, base2positions
    # mv only if finished
    os.replace(ofn+".fastq.gz_", ofn+".fastq.gz")
    return fn, basecount, mods2count, alphabet, symbol2modbase, canonical2mods, base2positions

def get_output_fnames(fnames, indir, outdir):
    """Return FastQ directory name and output names for FastQ files"""
    ofnames = []
    fqdir = os.path.join(outdir, "reads", os.path.basename(indir.rstrip(os.path.sep)))
    for fn in fnames:
        ofn = os.path.join(fqdir, fn[len(indir):].lstrip(os.path.sep))
        ofnames.append(ofn)
        # create out directories for FastQ files
        if not os.path.isdir(os.path.dirname(ofn)): os.makedirs(os.path.dirname(ofn))
    return fqdir, ofnames

def mod_encode(outdir, indirs, threads, basecall_group="", MaxModsPerBase=3, recursive=False):
    """Convert basecalled Fast5 into FastQ with base modification probabilities
    encoded as FastQ qualities.
    """
    fastq_dirs = []
    MaxPhredProb = get_MaxPhredProb(MaxModsPerBase)
    logger("Encoding modification info from %s directories...\n"%len(indirs))
    for indir in indirs:
        if recursive:
            fnames = sorted(map(str, Path(indir).rglob('*.fast5')))
        else:
            fnames = sorted(map(str, Path(indir).glob('*.fast5')))
        logger(" %s with %s Fast5 file(s)...\n"%(indir, len(fnames)))
        # no need to have more threads than input directories ;) 
        if threads > len(fnames):
            threads = len(fnames)
        # define imap, either pool of processes or map
        if threads>1:
            p = Pool(threads, maxtasksperchild=1)
            imap = p.imap_unordered 
        else:
            imap = map
        # load current data - this may cause problems if recursion was done...
        data = load_info(indir, recursive)
        # exit if not fnames nor modPhred.pkl
        if not fnames and not data:
            warning("[mod_encode][WARNING] No Fast5 files and no previously process data in %s\n"%indir)
            sys.exit(1)
            continue
        # get already processed files
        fast5 = {}
        if data: 
            # start from the beginning if different MaxModsPerBase used
            if data["MaxModsPerBase"] != MaxModsPerBase:
                info = "[mod_encode][WARNING] Previously you used --MaxModsPerBase %s, while now %s. Recomputing FastQ files..."
                warning(info%(data["MaxModsPerBase"], MaxModsPerBase))
            else:
                fast5 = data["fast5"]
                logger("  %s were processed earlier."%len(fast5), add_memory=0)
        # process files if not already processed (FastQ not present or FastQ is older than Fast5)
        fqdir, ofnames = get_output_fnames(fnames, indir, outdir)
        fastq_dirs.append(fqdir)
        args = [(fn, ofn, basecall_group, MaxPhredProb) for fn, ofn in zip(fnames, ofnames)
                if not os.path.isfile(ofn+".fastq.gz")]
        parser = imap(worker, args)
        for ii, (fn, basecount, mods2count, alphabet, symbol2modbase, canonical2mods, base2positions) in enumerate(parser, 1):
            # skip files without bases
            if not basecount: continue
            sys.stderr.write(" %s / %s  %s with %s bases. Detected mods: %s   \r"%(ii, len(args), os.path.basename(fn), basecount, str(mods2count)))
            data = load_info(fqdir, recursive)
            # store data
            if data:
                # either add new Fast5
                fast5 = data["fast5"]
                fast5[fn] = basecount
                fast5mod = data["fast5mod"]
                fast5mod[fn] = mods2count
            else:
                # or start pickle from scratch if doesn't exists
                fast5 = {fn: basecount}
                fast5mod = {fn: mods2count}
                moltype, bases = get_moltype(alphabet)
                nmods = sum(len(v) for k, v in canonical2mods.items())
                info = "  %s alphabet with %s modification(s) %s. symbol2modbase: %s"
                logger(info%(moltype, nmods, str(canonical2mods), str(symbol2modbase)), add_memory=0)
                # make sure --MaxModsPerBase is sufficiently large
                maxnmodsperbase = max(len(v) for k, v in canonical2mods.items())
                if maxnmodsperbase>MaxModsPerBase:
                    info = "[mod_encode][WARNING] Too many modifications per base (%s). \nPlease restart with --MaxModsPerBase %s or larger!"
                    warning(info%(maxnmodsperbase, maxnmodsperbase))
            # this keeps info on completed Fast5>FastQ this way
            dump_info(fqdir, alphabet, symbol2modbase, canonical2mods, base2positions, fast5, fast5mod,
                      MaxModsPerBase, MaxPhredProb)
        # report total number of bases for project
        data = load_info(fqdir, recursive)
        symbol2modbase = data["symbol2modbase"]
        totbases = sum(v for k, v in data["fast5"].items())
        # get total number of modifications
        mods2count = {}
        if "fast5mod" in data:
            for fn in data["fast5mod"]:
                for k, v in data["fast5mod"][fn].items():
                    if k not in mods2count: mods2count[k] = v
                    else: mods2count[k] += v
        modcount = ", ".join(("{:,} {} [{:7,.3%}]".format(c, symbol2modbase[m], c/totbases)
                              for m, c in mods2count.items()))
        logger("  {:,} bases saved in FastQ, of those: {}   ".format(totbases, modcount), add_memory=0)
        # close pool
        if threads>1: p.terminate()
    return fastq_dirs

def warning(info):
    def count(i=0):
        i+=1
    logger(info, add_timestamp=0, add_memory=0)

def memory_usage(childrenmem=True, div=1024.):
    """Return memory usage in MB including children processes"""
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / div
    if childrenmem:
        mem += resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / div
    return mem
        
def logger(info, add_timestamp=1, add_memory=1, out=sys.stderr):
    """Report nicely formatted stream to stderr"""
    info = info.rstrip('\n')
    memory = timestamp = ""
    if add_timestamp:
        timestamp = "[%s]"%str(datetime.now()).split(".")[0] #"[%s]"%datetime.ctime(datetime.now())
    if add_memory:
        memory = " [mem: %5.0f MB]"%memory_usage()
    out.write("%s %s%s\n"%(timestamp, info, memory))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--indirs", nargs="+", help="input directory with Fast5 files")
    parser.add_argument("-o", "--outdir", default="modPhred", help="output directory [%(default)s]")
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("-t", "--threads", default=8, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("--basecall_group",  default="", help="basecall group to use from Fast5 file [last basecalling]")
    parser.add_argument("--MaxModsPerBase", default=MaxModsPerBase, type=int, help=argparse.SUPPRESS)

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    #sys.stderr.write("Processing %s directories...\n"%len(o.indirs))
    mod_encode(o.outdir, o.indirs, o.threads, o.basecall_group, o.MaxModsPerBase, o.recursive)

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

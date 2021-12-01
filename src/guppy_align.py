#!/usr/bin/env python3
desc="""Align FastQ files from input directories and store them in outdir/minimap2.

More info at: https://github.com/lpryszcz/modPhred

Dependencies: minimap2, samtools

TO DO:
- encode modification information in BAM files!
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 21/06/2019
"""

import glob, os, pickle, subprocess, sys
from datetime import datetime
from pathlib import Path
from guppy_encode import VERSION, is_rna, logger, load_info

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["./", "src"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])

def run_minimap2(ref, fastq, outfn, threads=1, storeQuals=False, spliced=1,
                 mem=1, tmpdir="/tmp", sensitive=False): 
    """Run minimap2 and sort bam on-the-fly"""
    mode = ["-axmap-ont", ]
    if spliced:
        mode = ["-axsplice", "-k13"] #-uf 
    # -k7 -w5 -m20 -A3 -B1 - very sensitive alignment
    args1 = ["minimap2", "-t%s"%threads, ref] + fastq + mode
    if sensitive:
        args1 += ["-k7", "-w5", "-m20", "-A3", "-B1"]
    # add query FastQ files
    #args1 += glob.glob(fastq)
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=open(outfn[:-4]+".log", "w"))
    # add baseq to SAM
    if storeQuals:
        args12 = ["sam_add_qual.py", "-f", ] + " ".join(fastq).replace(".fastm.gz", ".fastq.gz").split()
        proc12 = subprocess.Popen(args12, stdin=proc1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        proc12 = proc1
    # samtools sort process
    args2 = ["samtools", "sort", "-@%s"%threads, "-m%sG"%mem, "-o%s"%outfn, "-"] #"-T /tmp/%s"%fastq, 
    proc2 = subprocess.Popen(args2, stdin=proc12.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # wait to finish and close
    proc2.wait()
    proc2.terminate()

def mod_align(indirs, ref, outdir, threads, recursive=False, storeQuals=False, bamdir="minimap2"):
    """Align *.fq.gz files from input dirs and store sorted .bam in outdir/minimap2"""
    logger("Aligning FastQ files from %s directories...\n"%len(indirs))
    # prepare output directory
    if not os.path.isdir(os.path.join(outdir, bamdir)):
        os.makedirs(os.path.join(outdir, bamdir))
    # assert same setting has been used between different runs
    compare_info(indirs)
    # process indirs
    for indir in indirs:
        # run alignment if BAM doesn't exist already
        outfn = os.path.join(outdir, bamdir, os.path.basename(indir)+".bam")
        if os.path.isfile(outfn):
            logger(" %s already present"%outfn)
            continue
        logger("  > %s\n"%outfn)
        if recursive:
            fq = sorted(map(str, Path(indir).rglob('*.fastm.gz')))
        else:
            fq = sorted(map(str, Path(indir).glob('*.fastm.gz')))
        run_minimap2(ref, fq, outfn, threads, storeQuals, spliced=is_rna(indir))
        # update dump info
        dump_updated_info(indir, outdir, outfn)

def dump_updated_info(indir, outdir, bamfn, fn="modPhred.pkl"):
    """Load data from dump and store updated info.

    Check for version differences.
    """
    # load info from pickle
    newdata = load_info(indir)
    # store info
    data = load_info(outdir)
    if not data:
        data = newdata
        data["bam"] = []
    # add bam file
    data["bam"].append(bamfn)
    data["fast5"].update(newdata["fast5"])
    # dump using pickle
    fname = os.path.join(outdir, fn)
    with open(fname, "wb") as out:
        pickle.dump(data, out, protocol=pickle.HIGHEST_PROTOCOL)    
        
def compare_info(indirs, warnKeys=["version",],
                 errKeys=["MaxModsPerBase", "symbol2modbase", "alphabet"]):
    """Assert MaxModsPerBase, symbol2modbase and alphabet are the same for all runs.

    Warn if differences between software versions for runs. 
    """
    error = 0
    warns, errors = (), ()
    for indir in indirs:
        data = load_info(indir)
        _warns = tuple(data[x] if x in data else "missing" for x in warnKeys)
        _errors = tuple(data[x] if x in data else "missing" for x in errKeys)
        if warns:
            if warns!=_warns:
                sys.stderr.write("[WARNING] Difference at %s: %s vs %s\n"%(indir, warns, _warns))
        if errors:
            if errors!=_errors:
                sys.stderr.write("[ERROR] Difference at %s: %s vs %s\n"%(indir, errors, _errors))
                error = 1
        warns, errors = _warns, _errors
    if error:
        sys.exit(1)
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--indirs", nargs="+", help="input directories with FastQ files")
    parser.add_argument("-t", "--threads", default=8, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-f", "--fasta", required=1, type=argparse.FileType('r'), help="reference FASTA file")
    parser.add_argument("-o", "--outdir", default="modPhred", help="output directory [%(default)s]")
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("-s", "--storeQuals", action='store_true', help="store base qualities in BAM [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    # we need name, not file object here
    o.fasta = o.fasta.name

    # get alignments
    mod_align(o.indirs, o.fasta, o.outdir, o.threads, o.recursive, o.storeQuals)
    
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

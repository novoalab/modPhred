#!/usr/bin/env python3
desc="""Convert basecalled Fast5 with modifications annotated by guppy v3.1.5+
to FastQ with modification probabilities encoded as FastQ qualities.

More info at: https://github.com/lpryszcz/modPhred

Dependencies: h5py

TO DO:
- drop low quality reads
- ignore low quality bases
- store FastQ only if needed (~33% increase in speed)
- merge guppy_encode & guppy_encode_live ie replacing minimap2 with mappy
- add polyA tail estimation
- SAM methylation encoding
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 20/06/2019
"""

import ast, glob, mappy, os, pysam, sys
import numpy as np
from collections import OrderedDict
from datetime import datetime
from multiprocessing import Pool # use pebble instead
from pathlib import Path
from basecall import init_args, start_guppy_server, get_basecall_client, basecall_and_align
from common import *

def get_MaxPhredProb(canonical2mods, QUALS=QUALS, min_mods_per_base=3):
    MaxModsPerBase = max(min_mods_per_base, max(map(len, canonical2mods.values())))
    MaxPhredProb = int(len(QUALS)/MaxModsPerBase)
    return MaxPhredProb

def get_phredmodprobs(seq, modbaseprobNorm, mods2count, base2positions,
                      canonical2mods, MaxPhredProb, a):
    """Return PHRED scaled probability of given base being modified 
    and mods2count dictionary
    """
    # prepare phred str as numpy array
    phredmodprobs = np.empty(len(seq), dtype='|S1')
    phredmodprobs[:] = QUALS[0]
    # get positions of each base that can be modified
    seq = np.array(list(seq)) # ~10x faster than iterating through seq and comparing to b
    base2idx = {b: np.argwhere(seq==b)[:, 0] # 25% faster than .flatten()
                for b in base2positions if len(base2positions[b])>1}
    # calculate PHRED scaled modification probability for each base that can be modified
    for b, ii in base2idx.items():
        # only calculate indices max if more than 1 modification per base
        if len(base2positions[b])>2:
            s, e = base2positions[b][1], base2positions[b][-1]+1
            probs = modbaseprobNorm[ii, s:e]
            indices = np.argmax(probs, axis=1)
            probs = np.max(probs, axis=1)
        else:
            probs = modbaseprobNorm[ii, base2positions[b][1]] 
            indices = np.zeros(probs.shape[0], dtype='int')
        # update PHRED - this takes ~46% of the function time - can we do it faster?
        phredmodprobs[ii] = QUALS[probs+indices*MaxPhredProb]
        # calculate stats
        # don't count supplementary alignments
        if a.is_secondary or a.is_supplementary: continue        
        for idx in range(len(base2positions[b])-1):
            mods2count[canonical2mods[b][idx]] += np.sum(probs[indices==idx]>=0.5*MaxPhredProb)
    return phredmodprobs.tobytes().decode(), mods2count

def get_sam_header(ref, coord_sorted=False, version='1.6'):
    """Return pysam.AlignmentHeader object as OderedDict
    
    If coord_sorted==True, then coordinate sorting is info is added
    as the first element of the header.
    """
    faidx = pysam.FastaFile(ref)
    h = OrderedDict()
    # add version and sort info if needed
    if coord_sorted: h["HD"] = {'VN': version, 'SO': 'coordinate'}
    # add references
    h.update(pysam.AlignmentHeader.from_references(faidx.references,
                                                   faidx.lengths).as_dict())
    return h

def get_mod_data(bam):
    """Return modification data"""
    sam = pysam.AlignmentFile(bam)
    header_dict = sam.header.as_dict()
    fnames = []
    basecount = 0
    for ri, rg_dict in enumerate(header_dict["RG"], 0):
        fnames.append(rg_dict["ID"])
        info = ast.literal_eval(rg_dict["DS"])
        basecount += info["basecount"]
        if not ri:
            mods2count = info["mods2count"]
        else:
            for m, c in info["mods2count"].items():
                mods2count[m] += c
    return fnames, basecount, mods2count, info

def encode_mods_in_bam(args):
    """Basecall (or parse basecalled Fast5 file), align 
    and encode modification probability in BAM as base qualities.
    """
    fn, bam, ref, rna, conf, oq_tag = args
    mods2count, symbol2modbase = {}, {}
    if os.path.isfile(bam):
        # load info from BAM
        fnames, basecount, mods2count, md = get_mod_data(bam)
        alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(md['base_mod_alphabet'], md['base_mod_long_names'])
        return bam, basecount, mods2count, symbol2modbase
    # get sample name
    sample_name = os.path.basename(fn).split(".")[0]
    i = basecount = 0
    algs = []
    # generate BAM header
    header_dict = get_sam_header(ref, coord_sorted=True)
    header = pysam.AlignmentHeader.from_dict(header_dict)
    for a, sig, move, modbaseprob, md in basecall_and_align(fn, header, rna, conf):
        #if i>=100: break
        # get mod info
        if not algs:
            alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(md['base_mod_alphabet'], md['base_mod_long_names'])
            mods2count = {m: 0 for mods in canonical2mods.values() for m in mods}
            MaxPhredProb = get_MaxPhredProb(canonical2mods)
        # get mod probabilities and normalise to MaxPhredProb
        #print(modbaseprob); break
        modbaseprobNorm = (MaxPhredProb/MaxProb*modbaseprob).astype('uint8')
        # get read sequence
        seq = a.get_forward_sequence()
        # and base probabilities matching the read sequence
        if rna: modbaseprobNorm = modbaseprobNorm[::-1]
        '''# not needed since mappy2sam always returns soft-clipped algs
        # trim hard-clipped bases from modbaseprobs
        if a.is_reverse: 
            se = a.cigar[0][1] if a.cigar[0][0]==5 else None
            ss = a.cigar[-1][1] if a.cigar[-1][0]==5 else 0
        else:
            ss = a.cigar[0][1] if a.cigar[0][0]==5 else 0
            se = a.cigar[-1][1] if a.cigar[-1][0]==5 else None
        modbaseprobNorm = modbaseprobNorm[ss:se]
        '''
        # store original qualities
        if oq_tag: a.set_tag(oq_tag, a.qual)
        # store read group
        a.set_tag("RG", sample_name)
        # get modprobs as qualities
        phredmodprobs, mods2count = get_phredmodprobs(seq, modbaseprobNorm, mods2count, base2positions, canonical2mods, MaxPhredProb, a)
        # orient them accordingly to read orientation
        if a.is_reverse: phredmodprobs = phredmodprobs[::-1]
        # and store
        a.qual = phredmodprobs
        algs.append(a)
        # don't count supplementary alignments
        if a.is_secondary or a.is_supplementary: continue
        i += 1
        basecount += len(seq)
    # don't write BAM if nothing aligned
    if not i: return bam, basecount, mods2count, symbol2modbase
    # add basecount and modinfo to header
    mod_data = {"MaxPhredProb": MaxPhredProb, "basecount": basecount,
                "mods2count": mods2count, 
                "base_mod_alphabet": md['base_mod_alphabet'], 
                "base_mod_long_names": md['base_mod_long_names'], 
                }
    header_dict["RG"] = [{"ID": sample_name, "DS": str(mod_data)}]
    # sort algs
    algs = sorted(algs, key=lambda a: (a.reference_name, a.pos))
    # and write in BAM
    with pysam.AlignmentFile(bam, "wb", header=header_dict) as out:
        for a in algs: out.write(a)
    return bam, basecount, mods2count, symbol2modbase

def get_intermediate_dir_and_files(fnames, indir, outdir):
    """Return temp directory name and output names for intermediate BAM files"""
    bams = []
    tempdir = os.path.join(outdir, "minimap2", os.path.basename(indir.rstrip(os.path.sep)))
    for fn in fnames:
        bam = os.path.join(tempdir, fn[len(indir):].lstrip(os.path.sep)+".bam")
        bams.append(bam)
        # create out directories for FastQ files
        if not os.path.isdir(os.path.dirname(bam)): os.makedirs(os.path.dirname(bam))
    return tempdir, bams

def mod_encode(outdir, indirs, fasta, threads, rna, config, host, port,
               recursive, device, timeout, oq_tag):
    """Convert basecalled Fast5 into FastQ with base modification probabilities
    encoded as FastQ qualities.
    """
    logger("Encoding modification info from %s directories...\n"%len(indirs))
    # skip if BAM files exists - this can be done before anything else
    # start guppy server if needed
    conf = (config, host, port)
    if host:
        guppy_proc, host, port = start_guppy_server(host, config, port, device)
        # try connecting to guppy server first
        conf = (config, host, port)
        client, _get_read_data = get_basecall_client(*conf)
    else:
        guppy_proc = None
        logger("We'll use basecall information from Fast5 files...")
    # load reference for mappy
    logger("Loading reference index from %s..."%fasta)
    aligner = mappy.Aligner(fasta, preset="spliced" if rna else "map-ont", k=13)
    # start pool of workers
    # it's important to initialise the pool with aligner object as it can't be pickled
    p = Pool(threads, initializer=init_args, initargs=(aligner, )) #, maxtasksperchild=1)
    bams = []
    for indir in indirs:
        fnames = sorted(map(str, Path(indir).rglob('*.fast5') if recursive
                            else Path(indir).glob('*.fast5')))
        logger(" %s with %s Fast5 file(s)...\n"%(indir, len(fnames)))
        # process files if not already processed (FastQ not present or FastQ is older than Fast5)
        tempdir, _bams = get_intermediate_dir_and_files(fnames, indir, outdir)
        bam = tempdir + ".bam"
        bams.append(bam)
        if os.path.isfile(bam):
            # load info from BAM
            fnames, basecount, mods2count, md = get_mod_data(bam)#; print(md)
            alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(md['base_mod_alphabet'], md['base_mod_long_names'])
            modcount = ", ".join(("{:,} {} [{:7,.3%}]".format(c, symbol2modbase[m], c/basecount)
                                  for m, c in mods2count.items()))
            logger("  {:,} bases, of those: {}   ".format(basecount, modcount), add_memory=0)
            continue
        # exit if not fnames and not BAM created
        if not fnames:
            logger("[mod_encode][WARNING] No Fast5 files and no previously process data in %s\n"%indir)
            sys.exit(1)
        # fn, bam, ref, rna, conf, oq_tag
        args = [(fn, _bam, fasta, rna, conf, oq_tag) for fn, _bam in zip(fnames, _bams)]
        _bams, basecount, mods2count = [], 0, {}
        for ii, (_bam, _basecount, _mods2count, symbol2modbase) in enumerate(p.imap(encode_mods_in_bam, args), 1):
            # skip files without bases
            if not _basecount: continue
            # update info
            _bams.append(_bam)
            basecount += _basecount
            for m, c in _mods2count.items():
                if m not in mods2count: mods2count[m] = c
                else: mods2count[m] += c
            sys.stderr.write(" %s / %s %s bases. Detected mods: %s   \r"%(ii, len(args), basecount, str(mods2count)))
        # write sample stats
        modcount = ", ".join(("{:,} {} [{:7,.3%}]".format(c, symbol2modbase[m], c/basecount)
                              for m, c in mods2count.items()))
        logger("  {:,} bases, of those: {}   ".format(basecount, modcount), add_memory=0)
        # merge bam files
        pysam.merge("--write-index", "-p", "-@ %s"%threads, bam, *_bams)
        # and clean-up
        os.system("rm -r %s"%tempdir)
    # close pool
    p.close()
    # and guppy basecall server if it was started by this process
    if guppy_proc: guppy_proc.terminate()
    return bams

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
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-t", "--threads", default=2, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("--tag", default="", help="SAM tag to store original qualities ie. OQ [skipped]")
    guppy = parser.add_argument_group("Basecalling options")
    guppy.add_argument("-c", "--config", default="dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg", help="guppy model [%(default)s]")
    guppy.add_argument("--host", "--guppy_basecall_server", default="",
                        help="guppy server hostname or path to guppy_basecall_server binary [use basecall information from Fast5]")
    guppy.add_argument("--port", default=5555, type=int,
                        help="guppy server port (this is ignored if binary is provided) [%(default)s]")
    guppy.add_argument("--device", default="cuda:0", help="CUDA device to use [%(default)s]")
    guppy.add_argument("--timeout", default=10*60, help="timeout in seconds to process each Fast5 file [%(default)s]")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    mod_encode(o.outdir, o.indirs, o.fasta, o.threads, o.rna, o.config, o.host, o.port,
               o.recursive, o.device, o.timeout, o.tag)

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

#!/usr/bin/env python3

import gzip, os, pyguppyclient, pysam, re, sys, subprocess, tempfile, time
from ont_fast5_api.fast5_interface import get_fast5_file
from pyguppyclient import GuppyClientBase, yield_reads
from pyguppyclient.ipc import SimpleRequestType, SimpleReplyType
from common import logger
from mapping import mappy, mappy2sam
from collections import OrderedDict

def init_args(*args):
    """Share globals with pool of workers. 

    Note, this function has to be present in the actual file 
    that require given global!
    """
    global aligner
    aligner, = args

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

def basecall_and_align(fn, header, rna, conf, MD=False, only_signal=False):
    """Basecall reads from Fast5 file, encode modifications, align to reference and yield
    - read alignment (pysam.AlignedSegment) to the reference
    - raw read signal from Fast5 file
    - move table returned by basecaller
    - modbaseprob table returned by basecaller

    All valid alignments are returned. 
    If MD==True, MD tag is additionally stored in alignment. 
    """
    # start client
    client = None
    if conf[1]: # host is defined
        client, _get_read_data = get_basecall_client(*conf)
        basecalled_reads = get_read_data_basecall(fn, client, _get_read_data)
    else: # read data from basecalled Fast5
        basecalled_reads = get_read_data(fn, only_signal=only_signal)
    # get index buffer from shared memory
    thr_buf = mappy.ThreadBuffer()
    # get basecalled reads from either guppy server or basecalled Fast5 file
    for read_id, seq, qseq, sig, move, modbaseprob, md in basecalled_reads:
        if rna: seq = seq.replace("U", "T")
        # iterate hits
        primary_count = 0
        for hit in aligner.map(seq, buf=thr_buf, MD=MD): 
            # catch secondary alignments
            if hit.is_primary: primary_count += 1
            a = mappy2sam(read_id, seq, qseq, hit, header, primary_count)
            yield a, sig, move, modbaseprob, md
    # disconnect guppy client
    if client: client.disconnect() # this is a bit dirty, but works fine

# Functions for basecalled Fast5
 
def __reduce_trace(trace):
    """Return trace summed for 4"""
    trace[:, :4] += trace[:, 4:] # this is way larger than raw_data
    trace = trace[:, :4] # this reduces the trace size (thus nearly all data) 2x
    return trace
        
def get_read_data(fn, only_signal=False):
    """Parse reads from basecalled Fast5 and yield: read_id
    qualities, raw signal, move and trace tables, and metadata. 
   
    If only_signal=True, it won't load move and trace tables, and metadata. 
    """
    read2data = {}
    f5file = get_fast5_file(fn, mode="r")
    for read in f5file.get_reads():
        read_id = read.read_id
        bcgrp = read.get_latest_analysis("Basecall_1D")
        fq = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Fastq")
        header, seq, _, qseq = fq.split("\n")[:4]
        sig = read.get_raw_data(scale=False) # True
        if only_signal:
            yield read_id, seq, qseq, sig, [], [], {}
        else:
            #trace = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Trace")
            modbaseprob = read.get_analysis_dataset(bcgrp, "BaseCalled_template/ModBaseProbs")
            move = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Move")
            bc_summary = read.get_analysis_attributes("%s/Summary/basecall_1d_template"%bcgrp)
            seggrp = read.get_latest_analysis("Segmentation")
            seg = read.get_analysis_attributes("%s/Summary/segmentation"%seggrp)
            mod_summary = read.get_analysis_attributes("%s/BaseCalled_template/ModBaseProbs"%bcgrp)
            channel_info = read.get_channel_info()
            md = {'trimmed_samples': seg['first_sample_template'],
                  'model_stride': bc_summary['block_stride'],
                  'mean_qscore': bc_summary['mean_qscore'],
                  'daq_offset': channel_info['offset'],
                  'daq_scaling': channel_info['range']/channel_info['digitisation'],
                  'base_mod_long_names': mod_summary['modified_base_long_names'].split(), #'6mA 5mC'
                  'base_mod_alphabet': mod_summary['output_alphabet'], #'AYCZGT'
                  }
            yield read_id, seq, qseq, sig, move, modbaseprob, md #__reduce_trace(trace)

# Functions for live basecalling

def get_read_data_basecall(fn, client, _get_read_data, batch_size=100):
    """Basecall reads from given Fast5 file and yield: read_id, sequence, 
    qualities, raw signal, move and trace tables. 

    The file is basecalled in batches of batch_size reads at once. 
    """
    # we need to count completed reads since v0.0.6 doesn't return raw data
    read2sig = {}
    # submit all reads
    for ri, read in enumerate(yield_reads(fn), 1):
        # zmq.error.Again: Resource temporarily unavailable
        try:
            client.pass_read(read)
        except:
            time.sleep(0.1)#; print('retrying pass_read()')
            client.pass_read(read)
        # store signal - v0.0.6 doesn't report raw data
        read2sig[read.read_id] = read.signal # read.daq_scaling*(read.signal+read.daq_offset)
        # grab basecalled reads in batches (initially there will be none)
        if ri%batch_size==0:
            # wait for the first batch to basecall
            if ri==len(read2sig): time.sleep(1)
            for read_id, seq, qseq, move, modbaseprob, md in _get_read_data(client):
                yield read_id, seq, qseq, read2sig.pop(read_id), move, modbaseprob, md
    # wait for the rest of the reads of basecalled reads
    while read2sig: 
        time.sleep(0.1)
        for read_id, seq, qseq, move, modbaseprob, md in _get_read_data(client):
            yield read_id, seq, qseq, read2sig.pop(read_id), move, modbaseprob, md
                
def start_guppy_server(host, config, port, device):
    """Start guppy_basecall_server and return its Popen object, host and port.

    This starts new basecall server only if host is a path to guppy binary. 
    Otherwise, it'll use host:port provided on the startup. 
    """
    def get_guppy_port(tmp, pat=re.compile(r'Starting server on port:\W+(\d+)')):
        """Return guppy port"""
        if os.path.isfile(tmp):
            for line in open(tmp):
                for m in pat.finditer(line): 
                    return int(m.groups()[0])
        return None
    
    proc = None
    if os.path.isfile(host):
        logger("Starting %s ..."%host)
        # get uniquely named tmp file
        tmp = os.path.join(tempfile.gettempdir(), next(tempfile._get_candidate_names()))
        # create log
        log = open(tmp+".log", "w")
        # start guppy process
        args = [host, "-l", tmp, "-c", config, "-x", device, "-p", "auto"]#; print(args)
        proc = subprocess.Popen(args, shell=False, stdout=log, stderr=log)
        while True:
            # get guppy port
            port = get_guppy_port(log.name)
            if port is not None: break
            # if no port, check if guppy init failed
            exitcode = proc.poll()
            if exitcode is not None:
                sys.stderr.write("There was some issue while starting guppy.\n")
                sys.stderr.write("Check logs for details: grep '' %s* \n"%tmp)
                # output guppy log to stderr
                log.close()
                for line in open(log.name): sys.stderr.write(line)
                sys.exit(exitcode)
            time.sleep(0.1)
        host = "localhost"
    return proc, host, port

def get_basecall_client(config, host, port):
    """Return client connection"""
    # define parameters for pyguppyclient v0.0.6 or newer
    ver = pyguppyclient.__version__
    kwargs = {} # no trace=True for v0.0.6, v0.0.7a1
    if ver=="0.0.6": _get_read_data = _get_read_data_v006
    elif ver=="0.0.7a1": _get_read_data = _get_read_data_v007a1
    elif ver in ("0.0.9", "0.1.0"):
        kwargs = {"trace": True}
        _get_read_data = _get_read_data_v009
    else: 
        sys.stderr.write("[ERROR] Unsupported pyguppy version: %s\n"%ver)
        sys.exit(1)
    # here due to v0.0.7a1 errors, we have to skip with and use connect after patching
    client = GuppyClientBase(config_name=config, host=host, port=port, **kwargs)
    # patch v0.0.7a1 that doesn't return trace and mod_probs
    if pyguppyclient.__version__=="0.0.7a1":
        client.pcl_client.set_params({'move_and_trace_enabled': True})
    client.connect()
    return client, _get_read_data

def _get_completed_reads(client, trace=True, state=False):
    """Get completed reads from pyguppyclient v0.0.6"""
    reads = []
    flag = (not trace) ^ state << 1
    res = client.send(SimpleRequestType.GET_FIRST_CALLED_BLOCK, data=flag)
    while res is not None:
        read, called = res
        while not called.complete:
            # zmq.error.Again: Resource temporarily unavailable
            try:
                _, block = client.send(SimpleRequestType.GET_NEXT_CALLED_BLOCK, data=read.read_tag)
            except:
                time.sleep(0.1)
                _, block = client.send(SimpleRequestType.GET_NEXT_CALLED_BLOCK, data=read.read_tag)                
            called += block
        # store read
        reads.append((read, called))
        res = client.send(SimpleRequestType.GET_FIRST_CALLED_BLOCK, data=flag)
    return reads

def _get_read_data_v006(client):
    """Yield read_id, seq, qseq, move and trace and metadata from baescalled reads"""
    for read, called in _get_completed_reads(client): 
        read_id = read.read_id
        seq, qseq = called.seq, called.qual
        move, trace = called.move, called.trace
        #trace = (trace * 255).astype("uint8") # called returns traces as 0-1 float
        modbaseprobs = (called.mod_probs*255).astype("uint8") # called returns modbaseprobs as 0-1 float
        md = {'trimmed_samples': called.trimmed_samples,
              'model_stride': called.model_stride,
              'mean_qscore' : called.qscore,
              'daq_offset': read.daq_offset, 
              'daq_scaling': read.daq_scaling,
              #'median': called.scaling["median"], 
              #'med_abs_dev': called.scaling["med_abs_dev"], 
              'base_mod_alphabet': called.mod_alphabet, # 'AYCZGT'
              'base_mod_long_names': called.mod_long_names, # ['6mA', '5mC']
              }
        # add scaling info
        md.update(called.scaling)
        yield read_id, seq, qseq, move, modbaseprobs, md

def __get_data_from_basecalled_read(read):
    """Return read_id, seq, qseq, move and trace and metadata from basecalled reads"""
    md, ds = read["metadata"], read["datasets"]
    read_id = md["read_id"]
    seq, qseq = ds["sequence"], ds["qstring"]
    move = ds["movement"]
    trace = ds["flipflop_trace"]
    modbaseprobs = ds["base_mod_probs"]
    return read_id, seq, qseq, move, modbaseprobs, md
        
def _get_read_data_v007a1(client):
    """Store read id, sequence and base_mod_probs. Compatible with v0.0.7a1."""
    for read in client.pcl_client.get_completed_reads()[0]:
        yield __get_data_from_basecalled_read(read)

def _get_read_data_v009(client):
    """Store read id, sequence and base_mod_probs. Compatible with v0.0.9."""
    for read in client.pcl_client.get_completed_reads():
        yield __get_data_from_basecalled_read(read)
        


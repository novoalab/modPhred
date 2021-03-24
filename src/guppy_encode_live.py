#!/usr/bin/env python3
desc="""Convert basecalled Fast5 with modifications annotated by guppy v3.1.5+
to FastQ with modification probabilities encoded as FastQ qualities.

More info at: https://github.com/lpryszcz/modPhred

Dependencies: h5py, pyguppyclient, running guppy_basecall_server

TO DO:
- catch exception of guppy client
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 17/02/2021
"""

from guppy_encode import *
import pyguppyclient, re, time, tempfile
from pyguppyclient import GuppyClientBase, yield_reads
from pyguppyclient.ipc import SimpleRequestType, SimpleReplyType

def get_completed_reads(client, trace=True, state=False):
    """Get completed reads from pyguppyclient v0.0.6"""
    reads = []
    flag = (not trace) ^ state << 1
    res = client.send(SimpleRequestType.GET_FIRST_CALLED_BLOCK, data=flag)
    while res is not None:
        read, called = res
        while not called.complete:
            _, block = client.send(SimpleRequestType.GET_NEXT_CALLED_BLOCK, data=read.read_tag)
            called += block
        # store read
        reads.append((read, called))
        res = client.send(SimpleRequestType.GET_FIRST_CALLED_BLOCK, data=flag)
    return reads

def _get_read_data_v006(client):
    """Store read id, sequence and base_mod_probs. Compatible with v0.0.6"""
    reads = []
    for read, called in get_completed_reads(client):
        reads.append((read.read_id, called.seq, called.mod_probs*255, 
                      " ".join(called.mod_long_names), called.mod_alphabet))
    return reads

def _get_read_data_v007a1(client):
    """Store read id, sequence and base_mod_probs. Compatible with v0.0.7a1. 
    """
    reads = [] # v0.0.7a1 returns reads, error_msg
    for read in client.pcl_client.get_completed_reads()[0]:
        md, ds = read["metadata"], read["datasets"]
        reads.append((md["read_id"], ds["sequence"], ds["base_mod_probs"],
                      md["base_mod_long_names"], md["base_mod_alphabet"]))
    return reads
            
def _get_read_data_v009(client):
    """Store read id, sequence and base_mod_probs. Compatible with v0.0.9. 
    """
    reads = [] # v0.0.9 returns reads
    for read in client.pcl_client.get_completed_reads():
        md, ds = read["metadata"], read["datasets"]
        reads.append((md["read_id"], ds["sequence"], ds["base_mod_probs"],
                      md["base_mod_long_names"], md["base_mod_alphabet"]))
    return reads
            
def get_basecalled_reads_data(fn, client, _get_read_data):
    """Return basecalled reads from given Fast5 file. 

    This implementation is >10x faster than using GuppyBasecallerClient.basecall()
    """
    reads = []
    # submit all reads
    for ri, read in enumerate(yield_reads(fn), 1):
        #if ri>100: break
        client.pass_read(read)
        # gradually grab basecalled reads (initially there will be none)
        if not ri%100: reads += _get_read_data(client)
    # wait for the rest of the reads
    while len(reads)<ri: 
        time.sleep(.1)
        reads += _get_read_data(client)
    return reads

def basecalling_worker(args):
    """Basecalling worker. 

    Basecalling is running as a separate worker process, so GPU is fully loaded. 
    Here read objects are small (10-100 Mb per Fast5), thus easy to pickle. 
    """
    fn, config, host, port = args
    # define parameters for pyguppyclient v0.0.6 or newer
    ver = pyguppyclient.__version__
    kwargs = {} # no trace=True for v0.0.6, v0.0.7a1
    if ver=="0.0.6": _get_read_data = _get_read_data_v006
    elif ver=="0.0.7a1": _get_read_data = _get_read_data_v007a1
    elif ver=="0.0.9":
        kwargs = {"trace": True}
        _get_read_data = _get_read_data_v009
    else: 
        sys.stderr.write("[ERROR] Unsupported pyguppy version: %s\n"%ver)
        sys.exit(1)
    # here due to v0.0.7a1 errors, we have to skip with and use connect after patching
    #with GuppyClientBase(config_name=config, host=host, port=port, **kwargs) as client:
    client = GuppyClientBase(config_name=config, host=host, port=port, **kwargs)
    # patch v0.0.7a1 that doesn't return trace and mod_probs
    if pyguppyclient.__version__=="0.0.7a1":
        client.pcl_client.set_params({'move_and_trace_enabled': True})
    client.connect()
    reads = get_basecalled_reads_data(fn, client, _get_read_data)
    client.disconnect() # this is a bit dirty, but works fine
    return fn, reads

def get_encoded_FastQ(reads, fn, MaxPhredProb):
    """Store modificoation probability in FastQ"""
    basecount = warns = 0
    data, rname = [], ""
    alphabet, symbol2modbase, canonical2mods, base2positions, mods2count = '', {}, {}, {}, {}
    # open out file with gzip compression
    outfn = fn+".fq.gz_"
    # The default mode is "rb", and the default compresslevel is 9.
    out = gzip.open(outfn, "wt")
    for ri, (read_name, seq, mod_probs, mods, output_alphabet) in enumerate(reads, 1):
        # prepare data storage if not already prepared
        if ri==1:
            alphabet, symbol2modbase, canonical2mods, base2positions = get_alphabet(output_alphabet, mods)
            rna = True if "U" in alphabet else False 
            mods2count = {m: 0 for mods in canonical2mods.values() for m in mods}
        # get mod probabilities and normalise to MaxPhredProb
        modbaseprobNorm = np.array(mod_probs / MaxProb * MaxPhredProb, dtype='uint8')
        # reverse only if RNA
        if rna: modbaseprobNorm = modbaseprobNorm[::-1]
        basecount += len(seq)
        # get modprobs as qualities
        phredmodprobs, mod2count = get_phredmodprobs(seq, modbaseprobNorm, mods2count, base2positions, canonical2mods, MaxPhredProb)
        out.write("@%s\n%s\n+\n%s\n"%(read_name, seq, phredmodprobs))
    # report number of bases
    if rname: name = "/".join(fn.split("/")[-4:])
    # mv only if finished
    os.replace(fn+".fq.gz_", fn+".fq.gz")
    return basecount, mods2count, alphabet, symbol2modbase, canonical2mods, base2positions

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

def mod_encode(indirs, threads, config, host, port, MaxModsPerBase=3,
               recursive=False, remove=False, device="cuda:0"):
    """Convert basecalled Fast5 into FastQ with base modification probabilities
    encoded as FastQ qualities.
    """
    # start guppy server if needed
    guppy_proc, host, port = start_guppy_server(host, config, port, device)
    fast5_dirs = set()
    MaxPhredProb = get_MaxPhredProb(MaxModsPerBase)
    logger("Encoding modification info from %s directories...\n"%len(indirs))
    for indir in indirs:
        if recursive:
            fnames = sorted(map(str, Path(indir).rglob('*.fast5')))
        else:
            fnames = sorted(map(str, Path(indir).glob('*.fast5')))
        # process & remove Fast5 files modified more than --remove minutes ago
        if remove:
            now = time.time()
            fnames = list(filter(lambda f: now-os.path.getmtime(f)>=remove*60, fnames))
        logger(" %s with %s Fast5 file(s)...\n"%(indir, len(fnames)))
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
        p = Pool(1) #, maxtasksperchild=1000)
        args = [(fn, config, host, port) for fn in fnames
                if not os.path.isfile(fn+".fq.gz")
                or os.path.getmtime(fn)>os.path.getmtime(fn+".fq.gz")]
        for ii, (fn, reads) in enumerate(p.imap(basecalling_worker, args), 1):
            # store modification probability in FastQ
            (basecount, mods2count, alphabet, symbol2modbase, canonical2mods,
             base2positions) = get_encoded_FastQ(reads, fn, MaxPhredProb)                
            # skip files without bases
            if not basecount: continue
            sys.stderr.write(" %s / %s  %s with %s bases. Detected mods: %s   \r"%(ii, len(fnames), os.path.basename(fn), basecount, str(mods2count)))
            data = load_info(indir, recursive)
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
            dump_info(indir, alphabet, symbol2modbase, canonical2mods, base2positions, fast5, fast5mod,
                      MaxModsPerBase, MaxPhredProb)
    # report total number of bases for project
    data = load_info(indir, recursive)
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
    p.close()
    # and guppy basecall server if it was started by this process
    if guppy_proc: guppy_proc.terminate()
    return fast5_dirs

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
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("-t", "--threads", default=6, type=int, help="number of cores to use [%(default)s]")
    #parser.add_argument("--basecall_group",  default="", help="basecall group to use from Fast5 file [last basecalling]")
    parser.add_argument("--MaxModsPerBase", default=MaxModsPerBase, type=int, help=argparse.SUPPRESS)
    #parser.add_argument("--remove", default=0, type=int, help="remove processed Fast5 files older than --remove minutes")
    parser.add_argument("-c", "--config", default="dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg", help="guppy model [%(default)s]")
    parser.add_argument("--host", "--guppy_basecall_server", default="localhost",
                        help="guppy server hostname or path to guppy_basecall_server binary [%(default)s]")
    parser.add_argument("-p", "--port", default=5555, type=int,
                        help="guppy server port (this is ignored if binary is provided) [%(default)s]")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    #sys.stderr.write("Processing %s directories...\n"%len(o.indirs))
    mod_encode(o.indirs, o.threads, o.config, o.host, o.port, o.MaxModsPerBase, o.recursive)

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

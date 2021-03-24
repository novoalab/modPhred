<img align="right" height="70" src="/doc/logo.png">

# modPhred

modPhred is a pipeline for **detection and annotation of DNA/RNA modifications from raw ONT data**.
To run it, you will need Fast5 files basecalled by guppy_basecaller v3.1.5+ with models trained to detect
DNA/RNA modifications of interest and FastA with reference sequences.
Some models are provided along with [guppy](https://community.nanoporetech.com/downloads).
You can train models to detect novel DNA/RNA modifications using [taiyaki](https://github.com/nanoporetech/taiyaki). 

<img align="right" width="350" src="/doc/pipeline.png">

The pipeline consists of four steps / modules ([more here](/doc#methods)):
- **1. modEncode**: Encoding modification probabilities in FastQ (mod_encode.py)
- **2. modAlign**: Build alignments keepind modification information in BAMs (get_align.py)
- **3. modReport** Extraction of RNA modification information (bedGraph) and QC reports (mod_report.py)
- **4. modAnalysis**: Plotting venn diagrams (mod_plot.py), co-occurrence of modifications(mod_correlation.py) and per-read clustering based on modification profiles (mod_cluster.py)

All these scripts can be run separately or as a pipeline by executing `modPhred/run`.  

## Running the pipeline
First, make sure you have [all dependencies installed in your system](/doc#installation). 

Then running modPhred is as easy as:
```bash
~/src/modPhred/run -f reference.fasta -o modPhred/projectName -i input_fast5_folder1 [input_fast5_folder2 ... input_fast5_folderN]
```

This will generate in the output directory `modPhred/projectName`:
- `.bed` - annotated positions with modifications as [bedMethyl-formatted](doc/#bedMethyl) files
  - mod.bed - combined report of positions with detected modifications in at least one of the samples
  - minimap2/*.bam.bed - modified sites reported for each run/sample separetely
- `mod.gz` - [internal format](doc#modgz) with all predicted positions with likely modifications
- `minimap2/*.bam` - alignments in BAM format and with encoded modifications.
Modification probabilities can be viewed directly in [IGV](doc#visualisation).
- `mod.gz.svg` - [QC plots](doc#qc-plots)

In addition, FastQ with encoded modifications as base qualities will be stored in as
`input_fast5_folder*/*.fast5.fq.gz` files. 

In order to see detailed description of program parameters, just execute it with `-h` / `--help`.
Make sure your Fast5 files are basecalled with
[guppy v3.1.5+ with models trained to detect modifications](/doc#how-to-check-if-my-fast5-files-are-basecalled-with-modifications).

For more usage examples, have a look in [test directory](test). 
You can find more information regarding internals of modPhred toolkit in the [documentation](/doc).

## Why using modPhred?
Cause why not! And seriously, it is:
- **free** ([MIT licensed](/LICENSE)) & **fast** ([5-6x faster than megalodon](/test#run-modphred))
- **easy-to-use** & **versatile**: will do all for you with just one command (or at least encoding of modifications in FastQ, alignments, detection of modified positions, QC & plotting...)
- **powerfull** & **space-optimised**: it stores the modification status inside FastQ/BAM
  - no external files/DBs needed
  - [you can visualise modification status of all bases of all reads in your favourite genome browser ie IGV](/doc#igv)
  - you can remove all Fast5 files to save disk space
- **visually attractive**: it produces nice plots (or at least not so bad so far... still working on it;) )

### What does modPhred stand for?
The tool stores base modification status (probability of base having various types of modifications)
encoded inside FastQ/BAM file instead of base quality (also called Phred scores),
thus **mod** (for modification) & **Phred** (for base quality).. Or something like that :P  
Initially, this tool was called **Pszczyna** ([pronounced ˈpʂt͡ʂɨna](https://forvo.com/word/pszczyna/)),
but since almost no one could pronounce or memorise it, 
we came up with much easier, yet so much less sexy name...
You can find more random stuff about Pszczyna [here](/doc#where-the-name-comes-from). 

## Getting help 
If you have any questions, issues or doubts, first please have a look in the [documentation](/doc#faq).
You'll find answers to most common questions and solutions to nearly all your problems.  
If above is not true for you, please [open new issue](/issue). 

## Citation 
[How to cite this tool?](/doc#citation)

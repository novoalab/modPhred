Program output
==============
If the program was execute command like

.. code-block:: bash

   ~/src/modPhred/run -f reference.fasta -o modPhred/projectName -i guppy_out/projectName/sample*/workspace

modPhred will generate in the output directory ``modPhred/projectName``:

* mod.gz - internal format with all predicted positions with likely modifications
* minimap2/*.bam - alignments in BAM format and with encoded modifications.
  One BAM file will be generated for every sample (directory) provided as input ``-i``.
  Modification probabilities can be :doc:`viewed directly in IGV <visualisation>`.
* .bed - annotated positions with modifications as bedMethyl-formatted files
  * mod.bed - combined report of positions with detected modifications in any of the samples
  * minimap2/*.bam.bed - modified sites reported for each run separetely
* mod.gz.svg - QC plots
  * and additional plots in plots/ directory

In addition, FastQ with encoded modifications as base qualities will be stored in
``guppy_out/projectName/sample*/workspace/*.fast5.fq.gz`` files.

Data formats
------------
While using modPhred it'll be good to familirised yourself with couple of data formats.

bedMethyl
^^^^^^^^^
modPhred reports information about modified positions in bedMethyl format. This is tab-delimited files, compatible with BED (positions are 0-based, half-open) with several additional fields:

#. Reference chromosome or scaffold
#. Start position in chromosome
#. End position in chromosome
#. Name of item - short modification name
#. Score - median modification probability scaled from 0-1000.
#. Strandedness, plus (+), minus (-), or unknown (.)
#. Start of where display should be thick (start codon) - same as 2.
#. End of where display should be thick (stop codon) - same as 3.
#. Color value (RGB) - different color to various modifications, although if more than 7 mods the colors may repeat. Color intensity depends on modification frequency (darker means modification is more frequent).
#. Coverage - number of reads at this position
#. Percentage of reads that show given modification at this position in the genome/transcriptome

For example, the output for ``NC_000913.3:1,061-1,253`` looks like this:

.. code-block:: bash

   NC_000913.3     1089    1090    5mC     903     +       1089    1090    0,139,0 454     55
   NC_000913.3     1091    1092    5mC     806     -       1091    1092    0,97,0  354     38
   NC_000913.3     1167    1168    6mA     839     +       1167    1168    0,0,155 468     61
   NC_000913.3     1168    1169    6mA     871     -       1168    1169    0,0,187 464     74
   NC_000913.3     1207    1208    5mC     806     +       1207    1208    0,108,0 431     42
   NC_000913.3     1209    1210    5mC     968     -       1209    1210    0,175,0 407     69

bedMethyl files can be visualised in many genome browsers ie IGV.

mod.gz
^^^^^^
This is internal, tab-delimited format that contain information about
all predicted positions with likely modifications.
Position are 1-based (similar to VCF format).

For example, the mod.gz for ``NC_000913.3:1,061-1,253`` region will look like that:

.. code-block:: bash

   ###
   # Welcome to modPhred (ver. 1.0b)!
   # 
   # Executed with: ~/src/modPhred/run -f ref/ECOLI.fa -o modPhred/PRJEB22772 -i guppy3.4.1/PRJEB22772/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145/workspace guppy3.4.1/PRJEB22772/MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711/workspace -t3
   #
   # For each bam file 4 values are stored for every position:
   # - depth of coverage (only positions with >=25 X in at least one sample are reported)
   # - accuracy of basecalling (fraction of reads having same base as reference, ignoring indels)
   # - frequency of modification (fraction of reads with modification above given threshold)
   # - median modification probability of modified bases (0-1 scaled). 
   #
   # If you have any questions, suggestions or want to report bugs,
   # please use https://github.com/lpryszcz/modPhred/issues.
   # 
   # Let's begin the fun-time with Nanopore modifications...
   ###
   chr     pos     ref_base        strand  mod     modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145.bam depth        modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145.bam basecall_accuracy    modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145.bam mod_frequency        modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145.bam median_mod_prob      modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711.bam depth        modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711.bam basecall_accuracy    modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711.bam mod_frequency        modPhred/PRJEB22772/minimap2/MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711.bam median_mod_prob
   NC_000913.3     244     C       -       5mC     444     0.910   0.014   0.806   120     0.958   0.050   0.581
   NC_000913.3     420     C       +       5mC     464     0.978   0.713   0.935   132     0.962   0.644   0.935
   NC_000913.3     422     C       -       5mC     351     0.604   0.328   0.806   103     0.621   0.369   0.839
   ... 
   NC_000913.3     1090    C       +       5mC     454     0.941   0.520   0.903   134     0.970   0.545   0.871
   NC_000913.3     1092    C       -       5mC     354     0.833   0.379   0.806   103     0.854   0.320   0.806
   NC_000913.3     1168    A       +       6mA     468     0.998   0.607   0.839   143     1.000   0.573   0.806
   NC_000913.3     1169    A       -       6mA     464     0.996   0.735   0.871   131     1.000   0.557   0.806
   NC_000913.3     1208    C       +       5mC     431     0.910   0.297   0.806   135     0.963   0.422   0.806
   NC_000913.3     1210    C       -       5mC     407     0.865   0.686   0.935   119     0.899   0.681   0.968



Why base Y is detected as modified, while model only reports modifications for X?
---------------------------------------------------------------------------------
Let's assume your model detects 5mC. Sometimes non-C reference bases may be detected as modified.
This may happend for several reasons:
* mis-alignment - apparent 5mC bases were incorrectly aligned to A, G or T reference
* mis-calling - apparent A, G or T bases were mispredicted as 5mC
* true biological variation - for example:
  * genotype of your sample may be different that this of your reference genome,
    thus true base will be C instead of A, G or T
  * heterozygous positions - a variant can have alternative allel being modified,
    thus 5mC may be true
  * variability in population - if you sequence pooled/mixed/tumor sample,
    some fraction of the cells may carry alternative alleles

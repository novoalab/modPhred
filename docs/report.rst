Detection of modifications
==========================
In order to detect and annotated modified positions, modPhred scans read alignments
and for every position of reference genome/transcriptome it calculates

* depth of coverage,
* basecall accuracy,
* modification frequency
* and median modification probability.

Those values are reported into ``mod.gz`` file.

By default, modPhred reports only sites with modifications fulfilling following criteria:

* minimum mapping quality of 15 - only reads with mapping quality of 15 or more are considered
* minimum sequencing depth of 25 - at least 25 reads for given positions
* at least 0.05 frequency of modification - at least 5% of reads being modified for this position
* at least 0.5 modification probability - only bases with modification probability of 50% or more are considered as truly modified.
  
The default values can be adjusted with several parameters:

.. code-block:: bash
		
   ./run -h
   ...
   -m --mapq               min mapping quality [15]
   -d --minDepth           min depth of coverage [25]
   --minModFreq            min modification frequency per position [0.05]
   --minModProb            min modification probability per base [0.5]
   ...


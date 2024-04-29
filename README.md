<img align="right" height="70" src="/docs/logo.png">

# modPhred

modPhred is a pipeline for detection, annotation and visualisation of DNA/RNA modifications
from raw ONT data.

The pipeline consists of four steps / modules:
1. **modEncode**: Encoding modification probabilities in FastQ (mod_encode.py)  
2. **modAlign**: Build alignments keeping modification information in BAMs (mod_align.py)  
3. **modReport** Extraction of RNA modification information (bedGraph) and QC reports (mod_report.py)  
4. **modAnalysis**: Plotting venn diagrams (mod_plot.py), co-occurrence of modifications (mod_correlation.py) and per-read clustering based on modification profiles (mod_cluster.py)  

For more information, please visit full documentation at https://modphred.readthedocs.io. 

## Does it work with any basecaller?

No, ModPhred works with Guppy flip-flip basecalling models (which can be custom trained, e.g. modification-aware basecalling models trained for example, with taiyaki). 
Unfortunately, ModPhred does not work with Dorado.

**As of March 2024, ONT has deprecated Guppy**. Consequently, as of March 2024, **ModPhred has also become deprecated.**
We will update this README if in the near future we upgrade ModPhred to support Dorado or alternative basecallers. 

## Citation 
If you find this work useful, please cite:

Pryszcz LP and Novoa EM (2022) ModPhred: an integrative toolkit for the analysis
and storage of nanopore sequencing DNA and RNA modification data.
[Bioinformatics, 38:257-260](http://dx.doi.org/10.1093/bioinformatics/btab539).


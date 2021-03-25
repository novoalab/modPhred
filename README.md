<img align="right" height="70" src="/docs/logo.png">

# modPhred

modPhred is a pipeline for detection, annotation and visualisation of DNA/RNA modifications
from raw ONT data. The pipeline consists of four steps / modules:

The pipeline consists of four steps / modules:
**1. modEncode**: Encoding modification probabilities in FastQ (mod_encode.py)
**2. modAlign**: Build alignments keepind modification information in BAMs (mod_align.py)
**3. modReport** Extraction of RNA modification information (bedGraph) and QC reports (mod_report.py)
**4. modAnalysis**: Plotting venn diagrams (mod_plot.py), co-occurrence of modifications(mod_correlation.py) and per-read clustering based on modification profiles (mod_cluster.py)

For more information, please visit full documentation at https://modphred.readthedocs.io. 

## Citation 
Manuscript is in preparation. Till preprint is ready, please cite this repository.
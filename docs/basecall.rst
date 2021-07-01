Basecalling
===========

Basecalling is performed using guppy basecaller.
You can choose to either: 

- basecall Fast5 before (guppy ver. >3.2 is supported)
- or run live basecalling (guppy ver. >3.6 is supported)

In the live mode, basecalling is performed with default settings. 
By default, CpG/dam/dcm model (``dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg``) is used.
This can be changed using ``-c / --config`` parameter. 
  
Note, if you choose to run live basecalling, you'll need to 
`install matching version of pyguppyclient <install.html#which-pyguppyclient-version-should-i-install>`_.

Modification-aware models
-------------------------
Modification aware-model should be used for basecalling. 
At this point, there are no publically available models for RNA modifications. 
The RNA model provided in this repository has been trained 
with *in vitro* transcribed molecules 
in which either all bases are modified or all are unmodified,
and because of that it isn't applicable to biological samples in which only some bases
are expected to be modified. 
You can find more details in the supplementary information of the manuscript.



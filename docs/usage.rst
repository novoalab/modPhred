Running the pipeline
====================

ModPhred can process DNA and RNA datasets.
The only difference between the two is enabling splice-aware alignments for RNA.
The type of dataset is dectected automatically from provided guppy config file ``-c / --config``.

The only required inputs are:

* reference FastA sequence
* path(s) containing Fast5 files

You can provide multiple input Fast5 folders -
those will be treated as separate runs/samples.
If you wish to treat multiple runs as one sample, place your Fast5 files in 1 folder.

ModPhred can be executed in three modes:

* local on-the-fly basecalling
* remote on-the-fly basecalling
* without basecalling (assuming the Fast5 files were basecalled before)

We strongly recommend to use on-the-fly basecalling because
it's a few times faster than running basecalling and modPhred separately.
Basecalling is **the slowest step** of entire process
and it produces **huge intermediate files**. 

In order to see detailed description of program parameters,
just execute it with ``-h / --help``.

Local on-the-fly basecalling
----------------------------
Here, we assume, that you have guppy already installed in you system. ModPhred will start guppy_basecall_server in the background and stop it when it isn't needed anymore.

All you need to do is to provide path to guppy_basecall_server (--host)

.. code-block:: bash

   ~/src/modPhred/run --host ~/src/ont-guppy_4.0.15/bin/guppy_basecall_server -f reference.fasta -o modPhred/projectName -i input_fast5_folder1 [input_fast5_folder2 ... input_fast5_folderN]

Alternatively, if guppy_basecall_server is already running in your machine, you can provide just its port using --host localhost --port.

Remote on-the-fly basecalling
-----------------------------
Here, we assume the guppy_basecall_server is already running in the remote machine. All you need to do is to provide IP address --host and port --port

.. code-block:: bash

   ~/src/modPhred/run --host 172.21.11.186 --port 5555 -f reference.fasta -o modPhred/projectName -i input_fast5_folder1 [input_fast5_folder2 ... input_fast5_folderN]

Without basecalling
-------------------
Make sure your Fast5 files are basecalled with guppy v3.1.5+ with models trained to detect modifications.

Running modPhred pipeline is as easy as:

.. code-block:: bash

   ~/src/modPhred/run -f reference.fasta -o modPhred/projectName -i input_fast5_folder1 [input_fast5_folder2 ... input_fast5_folderN]

For more usage examples, please have a look in :doc:`test dataset <test>`.

Processing (very) large datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are several ways of speeding up entire analysis for very large datasets.

* modEncode: process each sample or (or even subsets of each run) separately using guppy_encove_live.py. Ideally, each subset will be processed on dedicated GPU (local or remote). Here, providing more than 6 cores per job brings no improvement, since modEncode is primarily GPU-bound.
* modAlign: no much can be done, since every sample has to produce one BAM file.
* Beside, modAlign is by far the fastest step.
* modReport: process each chromsome (or even subsets of chromosome) as separate job. Make sure to provide as many cores as possible to each job.


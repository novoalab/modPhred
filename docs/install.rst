Installation
============

modPhred is written in Python3 and should work in most UNIX systems.

Make sure you install all programs listed below, before runnning the pipeline.

* first the repository

.. code-block:: bash
		
   mkdir ~/src
   cd src
   git clone https://github.com/novoalab/modPhred

* from `conda <https://bioconda.github.io/user/install.html#install-conda>`_ (use miniconda3!)

.. code-block:: bash

   conda install minimap2 samtools hdf5 wget

* from `pip <https://pypi.org/project/pip/>`_

.. code-block:: bash

   pip install h5py matplotlib pysam pandas seaborn

* guppy_basecaller has to be obtained from `Nanopore Tech. Software page <https://community.nanoporetech.com/downloads>`_
Alternatively, you can try `this for GPU <https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_4.0.15_linux64.tar.gz>`_
or `this for CPU <https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_4.0.15_linux64.tar.gz>`_ version.
For GPU basecalling to work, you'll need to install CUDA with NVIDIA drivers.
Check `my blog for instructions for Ubuntu 18.04 <https://medium.com/@lpryszcz/containers-with-cuda-support-5467f393649f>`_
or `NVIDIA CUDA website for other systems <https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html>`_.

* pyguppyclient (this will work with guppy v4.0 - v4.3)

.. code-block:: bash

   pip install pyguppyclient==0.0.7a1

Once you have all dependencies installed,
we recommend to try running it with :doc:`test dataset <test>`.
It'll be much easier to troubleshoot all potential issues this way. 
   
Which pyguppyclient version should I install?
---------------------------------------------

Guppy API was changing across version and unfortunately the newer versions are not back-compatible.
Therefore, you'll need to install pyguppyclient version matching the version of guppy basecaller.

=============== ===============
 Guppy version   pyguppyclient
                 version
=============== ===============
 >= 4.4 	 0.0.9                 
 >= 4.0 & <4.4   0.0.7a1         
 >= 3.4 & <4.0   0.0.6           
 < 3.4 	         not supported!        
=============== ===============


For example, if you intend to use guppy 4.0.15, you'll need to install pyguppyclient v0.0.7a1 as follows:

.. code-block:: bash

   pip install pyguppyclient==0.0.7a1

Note, only one version of pyguppyclient can be installed in your system. If you wish to use more than one version, you can install them using virtual environments as follows:

.. code-block:: bash

   python3 -m venv ~/src/venv/pyguppyclient006
   source ~/src/venv/pyguppyclient006/bin/activate
   pip install pyguppyclient==0.0.6 pysam pandas seaborn

And whenever you wish to switch to this version, just execute:

.. code-block:: bash

   source ~/src/venv/pyguppyclient006/bin/activate

Once you are finish with computation eihert close the terminal window
or execute ``deactivate``.

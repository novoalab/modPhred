# Set the base image to guppy image
FROM genomicpariscentre/guppy-gpu:3.6.1

# File Author / Maintainer
MAINTAINER Leszek Pryszcz <lpryszcz@crg.es>

ARG PACKAGE_VERSION=1.0b
ARG BUILD_PACKAGES="wget apt-transport-https"
ARG DEBIAN_FRONTEND=interactive
ARG PYVER=3.7

# miniconda
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    sh Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3 && \
    rm Miniforge3-Linux-x86_64.sh && \
    echo ". /opt/miniforge3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo ". /opt/miniforge3/etc/profile.d/mamba.sh" >> ~/.bashrc && \
    . ~/.bashrc && \
    echo "mamba activate python$PYVER" >> ~/.bashrc && \
    mamba activate && mamba create -y -n "python$PYVER" python=$PYVER && mamba activate python$PYVER && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# dependencies
RUN . ~/.bashrc && \
    mamba install -y minimap2 samtools hdf5 pysam && \
    pip install h5py matplotlib pandas seaborn scikit-learn

# version-specific dependencies
RUN . ~/.bashrc && \
    pip install pyguppyclient==0.0.6

# tool
RUN wget -q https://github.com/novoalab/modPhred/archive/refs/heads/main.tar.gz && \
    tar xpfz main.tar.gz && rm main.tar.gz && mv modPhred-main /opt/modPhred

# activate conda environment
RUN ln -s /opt/modPhred /opt/app
WORKDIR /opt/app
ENV PATH /opt/miniforge3/envs/python3.7/bin:$PATH
RUN ./run -h

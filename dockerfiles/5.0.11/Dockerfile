# Set the base image to guppy image
FROM genomicpariscentre/guppy-gpu:5.0.11

# File Author / Maintainer
MAINTAINER Leszek Pryszcz <lpryszcz@crg.es>

ARG PACKAGE_VERSION=1.0b
ARG BUILD_PACKAGES="wget apt-transport-https"
ARG DEBIAN_FRONTEND=interactive

# miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    sh Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    echo ". /opt/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate" >> ~/.bashrc && \
    . ~/.bashrc && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# dependencies
RUN . ~/.bashrc && \
    conda install -y minimap2 samtools hdf5 && \
    pip install h5py matplotlib pysam pandas seaborn sklearn

# version-specific dependencies
RUN . ~/.bashrc && \
    apt update && apt install build-essential -y && \
    pip install cython pyguppyclient==0.1.0

# tool
RUN wget -q  https://github.com/novoalab/modPhred/archive/refs/heads/main.tar.gz && \
    tar xpfz main.tar.gz && rm main.tar.gz && mv modPhred-main /opt/modPhred

# clean-up
RUN apt purge build-essential -y && \
    apt autoremove --purge --yes && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# activate conda environment
ENV PATH /opt/miniconda3/bin:$PATH
RUN /bin/bash -c "source activate"
# Set the base image to guppy image
FROM genomicpariscentre/guppy-gpu:3.6.1

# File Author / Maintainer
MAINTAINER Leszek Pryszcz <lpryszcz@crg.es>

ARG PACKAGE_VERSION=3.6.1
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
    pip install h5py matplotlib pysam pandas seaborn

# version-specific dependencies
RUN . ~/.bashrc && \
    pip install pyguppyclient==0.0.6

# tool
RUN wget -q https://github.com/novoalab/modPhred/archive/refs/heads/main.tar.gz && \
    tar xpfz main.tar.gz && rm main.tar.gz && mv modPhred-main /opt/modPhred

# activate conda environment
ENV PATH /opt/miniconda3/bin:$PATH
RUN /bin/bash -c "source activate"
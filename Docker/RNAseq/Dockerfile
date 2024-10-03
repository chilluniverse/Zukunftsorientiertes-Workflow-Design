# Basis-Image: Verwende ein Ubuntu-Image
FROM continuumio/anaconda3:2024.06-1

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    zip \
    unzip \
    ca-certificates \
    software-properties-common \
    apt-transport-https \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libffi-dev \
    libreadline-dev \
    libtiff5-dev \
    libjpeg-dev \
    libpng-dev \
    libx11-dev \
    libxt-dev \
    libopenblas-dev \
    perl \
    make \
    bzip2 \
    gnuplot \
    gawk \
    gfortran \
    libhdf5-dev \
    python3-dev \
    python3-pip \
    r-base \
    && apt-get clean

RUN conda update -n base -c defaults conda
RUN conda config --add channels anaconda
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# Installation von Python-Paketen
RUN conda install -y conda-forge::pandas=2.2.3

# Installation von R-Paketen, inkl. ggplot2 Version 3.5.1
RUN conda install -y conda-forge::r-ggplot2=3.5.1
RUN conda install -y conda-forge::r-biocmanager=1.30.25
RUN conda install -y bioconda::bioconductor-deseq2=1.42.0
RUN conda install -y bioconda::bioconductor-coseq=1.26.0

# Cleanup unnötiger Dateien
RUN conda clean --all -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Standard-Arbeitsverzeichnis
WORKDIR /data

# Kommando für den Start
CMD ["/bin/bash"]

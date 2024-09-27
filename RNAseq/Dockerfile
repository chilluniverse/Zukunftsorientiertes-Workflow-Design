# Basis-Image: Verwende ein Ubuntu-Image
FROM ubuntu:20.04

# Maintainer-Information
LABEL maintainer="dein_name@beispiel.com"

# Umgebungsvariablen für non-interactive Installation
ENV DEBIAN_FRONTEND=noninteractive

# Update und Installation der Abhängigkeiten
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
    gfortran \
    libhdf5-dev \
    python3-dev \
    python3-pip \
    r-base \
    && apt-get clean

# Installation von Bowtie2 und Samtools
RUN apt-get install -y bowtie2 samtools

# Installation von Python-Paketen
RUN pip3 install pandas

# Installation von R-Paketen
RUN R -e "install.packages(c('BiocManager', 'ggplot2', 'Cairo'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('DESeq2', 'coseq'))"

# Cleanup unnötiger Dateien
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Standard-Arbeitsverzeichnis
WORKDIR /data

# Kommando für den Start
CMD ["/bin/bash"]

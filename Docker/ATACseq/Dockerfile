# Basis-Image auswählen (hier ein Ubuntu-Image)
FROM ubuntu:20.04

# Umgebungsvariablen setzen
ENV DEBIAN_FRONTEND=noninteractive

# Notwendige Pakete installieren
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    zip \
    make \
    gcc \
    perl \
    python3 \
    python3-pip \
    libgd-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libexpat1-dev \
    git \
    && apt-get clean

# Arbeitsverzeichnis setzen
WORKDIR /data

# HOMER herunterladen und entpacken
RUN wget http://homer.ucsd.edu/homer/configureHomer.pl \
    && perl configureHomer.pl -install

# HOMER in den Pfad aufnehmen
ENV PATH="$PATH:/data/bin"


# Das dm6-Paket (Drosophila melanogaster) installieren
RUN perl configureHomer.pl -install dm3
RUN perl configureHomer.pl -install dm6

# Standard-Kommando setzen
CMD ["/bin/bash"]

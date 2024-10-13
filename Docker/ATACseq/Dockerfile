# Import Base Image
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
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

# Set Work Directory
WORKDIR /data

# Download HOMER and install
RUN wget http://homer.ucsd.edu/homer/configureHomer.pl \
    && perl configureHomer.pl -install

# add HOMER to PATH
ENV PATH="$PATH:/data/bin"

# Install dm3 and dm6 for HOMER
RUN perl configureHomer.pl -install dm3
RUN perl configureHomer.pl -install dm6

CMD ["/bin/bash"]
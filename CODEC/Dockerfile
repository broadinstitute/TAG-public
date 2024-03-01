FROM --platform=linux/amd64 ubuntu:20.04

LABEL maintainer="lining@broadinstitute.org"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
    && apt-get install -y software-properties-common \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y python3.8 python3.8-dev python3.8-venv python3-pip \
    && pip3 install pandas argparse numpy pysam

RUN apt-get install -y r-base r-base-dev

RUN apt-get install -y \
    git \
    wget \
    bwa \
    libssl-dev \
    g++ \
    zlib1g-dev \
    autoconf \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    build-essential \
    software-properties-common 

RUN apt-get update -qq
RUN apt-get install -y openjdk-11-jdk



# Clone the CODECsuite repository
RUN git clone --recursive https://github.com/broadinstitute/CODECsuite.git /CODECsuite


RUN wget https://github.com/Kitware/CMake/releases/download/v3.28.0-rc3/cmake-3.28.0-rc3-linux-x86_64.tar.gz \
    && tar -xzvf cmake-3.28.0-rc3-linux-x86_64.tar.gz \
    && mv cmake-3.28.0-rc3-linux-x86_64 /opt/cmake-3.28 \
    && ln -s /opt/cmake-3.28/bin/cmake /usr/local/bin/cmake \
    && ln -s /opt/cmake-3.28/bin/ctest /usr/local/bin/ctest \
    && ln -s /opt/cmake-3.28/bin/cpack /usr/local/bin/cpack \
    && rm cmake-3.28.0-rc3-linux-x86_64.tar.gz

RUN cd CODECsuite && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

COPY dependencies/ /dependencies/
COPY reference_files/ /reference_files/

RUN cp /dependencies/samtools-1.9/samtools /usr/bin/
RUN chmod +x /CODECsuite/snakemake/script/agg_log.py
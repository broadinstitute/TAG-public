# Use Eclipse Temurin (the modern replacement for OpenJDK)
FROM eclipse-temurin:11-jre

# Install basics and Samtools
# This base image is built on Ubuntu, so apt-get works perfectly
RUN apt-get update && apt-get install -y \
    samtools \
    wget \
    tar \
    pv \
    && rm -rf /var/lib/apt/lists/*

# Install BBTools (which contains bbduk.sh)
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download -O bbmap.tar.gz \
    && tar -xvzf bbmap.tar.gz \
    && rm bbmap.tar.gz

# Add BBMap to the PATH
ENV PATH="/opt/bbmap:${PATH}"

# Verify installation
RUN samtools --version && bbduk.sh --version

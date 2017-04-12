#################################################################
# Dockerfile
#
# Software:         ucsc userApps
# Software Version: v325
# Description:      ucsc userApps image for SciDAP
# Website:          http://hgdownload.cse.ucsc.edu/admin/exe/userApps.v325.src.tgz, http://scidap.com/
# Provides:         bedGraphToBigWig
# Base Image:       scidap/scidap:v0.0.1
# Build Cmd:        docker build --rm -t scidap/ucsc-userapps:v325 .
# Pull Cmd:         docker pull scidap/ucsc-userapps:v325
# Run Cmd:          docker run --rm scidap/ucsc-userapps:v325 bash
#################################################################

### Base Image
FROM scidap/scidap:v0.0.1
MAINTAINER Andrey V Kartashov "porter@porter.st"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Install required packages (samtools)

RUN apt-get clean all &&\
    apt-get update &&\
    apt-get install -y --no-install-recommends \
    libssl-dev \
    libmysqlclient-dev \
    libpng12-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*

### Installing ucsc userApps

ENV VERSION 325
ENV NAME "ucsc-userapps"
ENV URL "http://hgdownload.cse.ucsc.edu/admin/exe/userApps.v${VERSION}.src.tgz"

RUN wget -q -O - $URL | tar -zxv --strip-components=2 && \
make -j 4 && \
rm ./bin/hg* && \
strip ./bin/*; true && \
cp ./bin/*Wig* /usr/local/bin/ && \
cp ./bin/bedTo* /usr/local/bin/ && \
cp ./bin/lift* /usr/local/bin/ && \
rm -rf ./*

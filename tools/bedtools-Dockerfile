#################################################################
# Dockerfile
#
# Software:         bedtools2
# Software Version: 2.25.0-pefs2
# Description:      bedtools2 image for SciDAP
# Website:          http://bedtools.readthedocs.org/, http://scidap.com/
# Provides:         bedtools2
# Base Image:       scidap/scidap:v0.0.1
# Build Cmd:        docker build --rm -t scidap/bedtools2:v2.25.0-pefs2 .
# Pull Cmd:         docker pull scidap/bedtools2:v2.25.0-pefs2
# Run Cmd:          docker run --rm scidap/bedtools2:v2.25.0-pefs2 bedtools2
#################################################################

### Base Image
FROM scidap/scidap:v0.0.1
MAINTAINER Andrey V Kartashov "porter@porter.st"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing bedtools2

ENV VERSION 2.25.0-pefs2
ENV NAME bedtools2
#ENV URL "https://github.com/arq5x/bedtools2/releases/download/v${VERSION}/bedtools-${VERSION}.tar.gz"
ENV URL "https://github.com/SciDAP/bedtools2/archive/v${VERSION}.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
cd ${NAME}-${VERSION} && \
make -j 4 && \
cd .. && \
cp ./${NAME}-${VERSION}/bin/bedtools /usr/local/bin/ && \
strip /usr/local/bin/*; true && \
rm -rf ./${NAME}-${VERSION}/

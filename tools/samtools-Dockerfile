#################################################################
# Dockerfile
#
# Software:         samtools
# Software Version: 1.2-242-4d56437
# Description:      samtools image for SciDAP
# Website:          https://samtools.github.io, http://scidap.com/
# Provides:         samtools/htslib/tabix/bgzip
# Base Image:       scidap/scidap:v0.0.1
# Build Cmd:        docker build --rm -t scidap/samtools:v1.2-242-4d56437 .
# Pull Cmd:         docker pull scidap/samtools:v1.2-242-4d56437
# Run Cmd:          docker run --rm scidap/samtools:v1.2-242-4d56437 samtools
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
    apt-get install -y  \
        libncurses5-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*

### Installing samtools/htslib/tabix/bgzip

ENV VERSIONH 1.2.1-254-6462e34
ENV NAMEH htslib
ENV URLH "https://github.com/samtools/htslib/archive/${VERSIONH}.tar.gz"
ENV SHA1H "6462e349d16e83db8647272e4b98d2a92992071f"

ENV VERSION 1.2-242-4d56437
ENV NAME "samtools"
ENV URL "https://github.com/samtools/samtools/archive/${VERSION}.tar.gz"
ENV SHA1 "4d56437320ad370eb0b48c204ccec7c73f653393"

RUN git clone https://github.com/samtools/htslib.git && \
cd ${NAMEH} && \
git reset --hard ${SHA1H} && \
make -j 4 && \
cd .. && \
cp ./${NAMEH}/tabix /usr/local/bin/ && \
cp ./${NAMEH}/bgzip /usr/local/bin/ && \
cp ./${NAMEH}/htsfile /usr/local/bin/ && \
#RUN wget -q -O - $URLH | tar -zxv && \
#cd ${NAMEH}-${VERSIONH} && \
#make -j 4 && \
#cd .. && \
#cp ./${NAMEH}-${VERSIONH}/tabix /usr/local/bin/ && \
#cp ./${NAMEH}-${VERSIONH}/bgzip /usr/local/bin/ && \
#cp ./${NAMEH}-${VERSIONH}/htsfile /usr/local/bin/ && \
strip /usr/local/bin/tabix; true && \
strip /usr/local/bin/bgzip; true && \
strip /usr/local/bin/htsfile; true && \
#ln -s ./${NAMEH}-${VERSIONH}/ ./${NAMEH} && \

git clone https://github.com/samtools/samtools.git && \
cd ${NAME} && \
git reset --hard ${SHA1} && \
make -j 4 && \
cp ./${NAME} /usr/local/bin/ && \
cd .. && \
strip /usr/local/bin/${NAME}; true && \
rm -rf ./${NAMEH}/ && \
rm -rf ./${NAME}/ && \
rm -rf ./${NAMEH}

#wget -q -O - $URL | tar -zxv && \
#cd ${NAME}-${VERSION} && \
#make -j 4 && \
#cd .. && \
#cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
#strip /usr/local/bin/${NAME}; true && \
#rm -rf ./${NAMEH}-${VERSIONH}/ && \
#rm -rf ./${NAME}-${VERSION}/

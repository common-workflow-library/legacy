class: DockerRequirement
dockerPull: scidap/alea:v1.2.2
#dockerImageId: scidap/alea:v1.2.2 #not yet ready
dockerFile: |
  #################################################################
  # Dockerfile
  #
  # Software:         Alea
  # Software Version: 1.2.2
  # Description:      Alea image for SciDAP
  # Website:          http://scidap.com/
  # Provides:         alea|bwa|Plink|vcftools|bedGraphToBigWig|shapeit2
  # Base Image:       scidap/scidap:v0.0.1
  # Build Cmd:        docker build --rm -t scidap/alea:v1.2.2 .
  # Pull Cmd:         docker pull scidap/alea:v1.2.2
  # Run Cmd:          docker run --rm scidap/alea:v1.2.2 alea
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


  ### Installing BWA

  ENV VERSION 0.7.12
  ENV NAME bwa
  ENV URL "https://github.com/lh3/bwa/archive/${VERSION}.tar.gz"

  RUN wget -q -O - $URL | tar -zxv && \
      cd ${NAME}-${VERSION} && \
      make -j 4 && \
      cd .. && \
      cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
      strip /usr/local/bin/${NAME}; true && \
      rm -rf ./${NAME}-${VERSION}/


  #  ### Installing bowtie
  #
  #  ENV VERSION 1.1.2
  #  ENV NAME bowtie
  #  ENV URL "https://github.com/BenLangmead/bowtie/archive/v${VERSION}.tar.gz"
  #
  #  RUN wget -q -O - $URL | tar -zxv && \
  #    cd ${NAME}-${VERSION} && \
  #    make -j 4 && \
  #    cd .. && \
  #    cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
  #    cp ./${NAME}-${VERSION}/${NAME}-* /usr/local/bin/ && \
  #    strip /usr/local/bin/*; true && \
  #    rm -rf ./${NAME}-${VERSION}/


  ### Installing samtools/htslib/tabix/bgzip

  ENV VERSIONH 1.2.1
  ENV NAMEH htslib
  ENV URLH "https://github.com/samtools/htslib/archive/${VERSIONH}.tar.gz"

  ENV VERSION "1.2"
  ENV NAME "samtools"
  ENV URL "https://github.com/samtools/samtools/archive/${VERSION}.tar.gz"

  RUN wget -q -O - $URLH | tar -zxv && \
      cd ${NAMEH}-${VERSIONH} && \
      make -j 4 && \
      cd .. && \
      cp ./${NAMEH}-${VERSIONH}/tabix /usr/local/bin/ && \
      cp ./${NAMEH}-${VERSIONH}/bgzip /usr/local/bin/ && \
      cp ./${NAMEH}-${VERSIONH}/htsfile /usr/local/bin/ && \
      strip /usr/local/bin/tabix; true && \
      strip /usr/local/bin/bgzip; true && \
      strip /usr/local/bin/htsfile; true && \
      ln -s ./${NAMEH}-${VERSIONH}/ ./${NAMEH} && \
      wget -q -O - $URL | tar -zxv && \
      cd ${NAME}-${VERSION} && \
      make -j 4 && \
      cd .. && \
      cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
      strip /usr/local/bin/${NAME}; true && \
      rm -rf ./${NAMEH}-${VERSIONH}/ && \
      rm -rf ./${NAME}-${VERSION}/


  ### Installing bedGraphToBigWig

  RUN  wget -q -O /usr/local/bin/bedGraphToBigWig  http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
  chmod a+x /usr/local/bin/bedGraphToBigWig


  ### Installing VCFTools

  ENV VERSION 0.1.14
  ENV NAME vcftools
  ENV URL "https://github.com/vcftools/vcftools/releases/download/v${VERSION}/${NAME}-${VERSION}.tar.gz"

  RUN wget -q -O - $URL | tar -zxv && \
      cd ${NAME}-${VERSION} && \
      ./configure --prefix=/usr/local && \
      make -j 4 install && \
      cd .. && \
      rm -rf ./${NAME}-${VERSION}/


  ### Installing SHAPEIT2

  RUN wget -q -O - https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz|tar -zxv -C /usr/local bin/shapeit && \
      chmod a+x /usr/local/bin/shapeit


  ### Installing plink

  RUN wget -q -O plink-1.07-x86_64.zip http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-x86_64.zip && \
      unzip plink-1.07-x86_64.zip && \
      cp plink-1.07-x86_64/* /usr/local/bin && \
      rm -rf plink-1.07-x86_64 && \
      rm -f plink-1.07-x86_64.zip


  ### Installing alea

  RUN wget -q -O - ftp://ftp.bcgsc.ca/supplementary/ALEA/files/alea.1.2.2.tar.gz | tar -zxv -C /usr/local --strip-components=1 && \
      cd /usr/local/bin/ && \
      printf '144c144\n<                 --output-fasta="$VAR_GENOME1_SNPS"\n---\n>                 --output-fasta="$VAR_GENOME2_SNPS"\n'| patch createGenome.sh && \
      sed -i.bak s/^AL_BWA_ALN_PARAMS/#AL_BWA_ALN_PARAMS/g alea.config && \
      sed -i.bak s/^AL_USE_CONCATENATED_GENOME/#AL_USE_CONCATENATED_GENOME/g alea.config && \
      rm -f alea.config.bak

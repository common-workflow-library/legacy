class: DockerRequirement
dockerPull: scidap/bcftools:v1.3
#dockerImageId: scidap/alea:v1.2.2 #not yet ready
dockerFile: |
  #################################################################
  # Dockerfile
  #
  # Software:         bcftools
  # Software Version: 1.3
  # Description:      bcftools image for SciDAP
  # Website:          http://samtools.github.io/bcftools/, http://scidap.com/
  # Provides:         bcftools/htslib/tabix/bgzip
  # Base Image:       scidap/scidap:v0.0.1
  # Build Cmd:        docker build --rm -t scidap/bcftools:v1.3 .
  # Pull Cmd:         docker pull scidap/bcftools:v1.3
  # Run Cmd:          docker run --rm scidap/bcftools:v1.3 bcftools
  #################################################################

  ### Base Image
  FROM scidap/scidap:v0.0.1
  MAINTAINER Andrey V Kartashov "porter@porter.st"
  ENV DEBIAN_FRONTEND noninteractive

  ################## BEGIN INSTALLATION ######################

  WORKDIR /tmp

  ### Installing bcftools/htslib/tabix/bgzip

  ENV NAMEH htslib
  ENV SHA1H "87141ea6ec0ee34db0748b2ae6a5840c456f5a02"

  ENV NAME bcftools
  ENV SHA1 "d388dd1aaef543a5ad1bf39454b657e8fe02fc80"

  RUN git clone https://github.com/samtools/htslib.git && \
      cd ${NAMEH} && \
      git reset --hard ${SHA1H} && \
      make -j 4 && \
      cd .. && \
      cp ./${NAMEH}/tabix /usr/local/bin/ && \
      cp ./${NAMEH}/bgzip /usr/local/bin/ && \
      cp ./${NAMEH}/htsfile /usr/local/bin/ && \
      strip /usr/local/bin/tabix; true && \
      strip /usr/local/bin/bgzip; true && \
      strip /usr/local/bin/htsfile; true && \

      git clone https://github.com/samtools/bcftools.git && \
      cd ${NAME} && \
      git reset --hard ${SHA1} && \
      make -j 4 && \
      cp ./${NAME} /usr/local/bin/ && \
      cp ./plugins/*.so /usr/local/bin/ && \
      cd .. && \
      strip /usr/local/bin/${NAME}; true && \
      rm -rf ./${NAMEH}/ && \
      rm -rf ./${NAME}/ && \
      rm -rf ./${NAMEH}

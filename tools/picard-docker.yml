class: DockerRequirement
dockerPull: scidap/picard:v1.141
#dockerImageId: scidap/picard:v1.141 #not yet ready
dockerFile: |
  #################################################################
  # Dockerfile
  #
  # Software:         picard
  # Software Version: 1.141
  # Description:      picard image for SciDAP
  # Website:          https://github.com/broadinstitute/picard/, http://scidap.com/
  # Provides:         picard

  # Base Image:       scidap/scidap:v0.0.1
  # Build Cmd:        docker build --rm -t scidap/picard:v1.141 .
  # Pull Cmd:         docker pull scidap/picard:v1.141
  # Run Cmd:          docker run --rm scidap/picard:v1.141 picard
  #################################################################

  ### Base Image
  FROM scidap/scidap:v0.0.1
  MAINTAINER Andrey V Kartashov "porter@porter.st"
  ENV DEBIAN_FRONTEND noninteractive

  ################## BEGIN INSTALLATION ######################

  WORKDIR /tmp

  ### INSTALL PICARD

  ENV VERSION "1.141"
  ENV NAME "picard-tools"
  ENV ZIP ${NAME}-${VERSION}.zip
  ENV URL https://github.com/broadinstitute/picard/releases/download/${VERSION}/${ZIP}

  RUN wget -q $URL -O ${ZIP} && \
      unzip $ZIP && \
      rm $ZIP && \
      cd ${NAME}-${VERSION} && \
      mv * /usr/local/bin && \
      cd .. && \
      bash -c 'echo -e "#!/bin/bash\njava -jar /usr/local/bin/picard.jar \$@" > /usr/local/bin/picard' && \
      chmod +x /usr/local/bin/picard && \
      rm -rf ${NAME}-${VERSION}



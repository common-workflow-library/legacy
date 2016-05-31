class: DockerRequirement
dockerPull:
dockerFile: |
  #################################################################
  # Dockerfile
  #
  # Software:         cnvkit
  # Software Version: 0.7.11
  # Description:      cnvkit docker image
  # Website:          https://github.com/etal/cnvkit
  # Provides:
  # Base Image:
  # Build Cmd:
  # Pull Cmd:
  # Run Cmd:
  #################################################################

  FROM python:2.7
  MAINTAINER Anton Khodak <anton.khodak@ukr.net>

  # Install cnvkit from pip
  RUN pip install cnvkit

  # Default command to execute at startup of the container
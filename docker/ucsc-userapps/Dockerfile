FROM ubuntu:15.04
MAINTAINER Andrey V Kartashov "porter@porter.st"
COPY ./ucsc-userapps/bin/*Wig* /usr/local/bin/
RUN apt-get update && apt-get install -y --no-install-recommends \
    libssl1.0.0 \
    && rm -rf /var/lib/apt/lists/*

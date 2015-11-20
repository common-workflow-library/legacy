#!/usr/bin/env bash

BASE_CHRS="\
chr1 \
chr2 \
chr3 \
chr4 \
chr5 \
chr6 \
chr7 \
chr8 \
chr9 \
chr10 \
chr11 \
chr12 \
chr13 \
chr14 \
chr15 \
chr16 \
chr17 \
chr18 \
chr19 \
chrX \
chrY \
chrM"
CHRS_TO_INDEX=$BASE_CHRS

UCSC_MM10_BASE=ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes

get() {
  file=$1
  if ! wget --version >/dev/null 2>/dev/null ; then
    if ! curl --version >/dev/null 2>/dev/null ; then
      echo "Please install wget or curl somewhere in your PATH"
      exit 1
    fi
    curl -o `basename $1` $1
    return $?
  else
    wget $1
    return $?
  fi
}

#INPUTS=
for c in $CHRS_TO_INDEX ; do
  if [ ! -f ${c}.fa ] ; then
    F=${c}.fa.gz
    get ${UCSC_MM10_BASE}/$F || (echo "Error getting $F" && exit 1)
    gunzip $F -c >>./mm10.fa || (echo "Error unzipping $F" && exit 1)

  fi
#  [ -n "$INPUTS" ] && INPUTS=$INPUTS,${c}.fa
#  [ -z "$INPUTS" ] && INPUTS=${c}.fa
done


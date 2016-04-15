#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: bedtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#i'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#sizeA'
    type:
      - 'null'
      - boolean
    description: "\t\tSort by feature size in ascending order.\n"
    inputBinding:
      position: 1
      prefix: '-sizeA'
  - id: '#sizeD'
    type:
      - 'null'
      - boolean
    description: "\t\tSort by feature size in descending order.\n"
    inputBinding:
      position: 1
      prefix: '-sizeD'
  - id: '#chrThenSizeA'
    type:
      - 'null'
      - boolean
    description: "\tSort by chrom (asc), then feature size (asc).\n"
    inputBinding:
      position: 1
      prefix: '-chrThenSizeA'
  - id: '#chrThenSizeD'
    type:
      - 'null'
      - boolean
    description: "\tSort by chrom (asc), then feature size (desc).\n"
    inputBinding:
      position: 1
      prefix: '-chrThenSizeD'
  - id: '#chrThenScoreA'
    type:
      - 'null'
      - boolean
    description: "\tSort by chrom (asc), then score (asc).\n"
    inputBinding:
      position: 1
      prefix: '-chrThenScoreA'
  - id: '#chrThenScoreD'
    type:
      - 'null'
      - boolean
    description: "\tSort by chrom (asc), then score (desc).\n"
    inputBinding:
      position: 1
      prefix: '-chrThenScoreD'
  - id: '#faidx'
    type:
      - 'null'
      - boolean
    description: "(names.txt)\tSort according to the chromosomes declared in \"names.txt\"\n"
    inputBinding:
      position: 1
      prefix: '-faidx'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: |
      Print the header from the A file prior to results.
    inputBinding:
      position: 1
      prefix: '-header'
outputs:
  - id: '#stdoutfile'
    type: File
    outputBinding:
      glob:
        engine: 'cwl:JsonPointer'
        script: /job/stdoutfile
stdout:
  engine: 'cwl:JsonPointer'
  script: /job/stdoutfile
baseCommand:
  - bedtools
  - sort
description: |
  Tool:    bedtools sort (aka sortBed)
  Version: v2.25.0
  Summary: Sorts a feature file in various and useful ways.

  Usage:   bedtools sort [OPTIONS] -i <bed/gff/vcf>

  Options: 
  	-sizeA			Sort by feature size in ascending order.
  	-sizeD			Sort by feature size in descending order.
  	-chrThenSizeA		Sort by chrom (asc), then feature size (asc).
  	-chrThenSizeD		Sort by chrom (asc), then feature size (desc).
  	-chrThenScoreA		Sort by chrom (asc), then score (asc).
  	-chrThenScoreD		Sort by chrom (asc), then score (desc).
  	-faidx (names.txt)	Sort according to the chromosomes declared in "names.txt"
  	-header	Print the header from the A file prior to results.


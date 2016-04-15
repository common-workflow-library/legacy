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
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Force strandedness.  That is, only merge features
      that are the same strand.
      - By default, merging is done without respect to strand.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#d'
    type:
      - 'null'
      - int
    description: |
      Maximum distance between features allowed for features
      to be merged.
      - Def. 0. That is, overlapping & book-ended features are merged.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-d'
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
  - cluster
description: |
  Tool:    bedtools cluster
  Version: v2.25.0
  Summary: Clusters overlapping/nearby BED/GFF/VCF intervals.

  Usage:   bedtools cluster [OPTIONS] -i <bed/gff/vcf>

  Options: 
  	-s	Force strandedness.  That is, only merge features
  		that are the same strand.
  		- By default, merging is done without respect to strand.

  	-d	Maximum distance between features allowed for features
  		to be merged.
  		- Def. 0. That is, overlapping & book-ended features are merged.
  		- (INTEGER)


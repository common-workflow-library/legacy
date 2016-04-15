#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: samtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#]'
    type: boolean
    description: ']'
    inputBinding:
      position: 4
  - id: '#in1bam'
    type: File
    description: in1.bam
    inputBinding:
      position: 2
  - id: '#b'
    type:
      - 'null'
      - boolean
    description: |
      <bed>            list of positions or regions
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#f'
    type:
      - 'null'
      - boolean
    description: |
      <list>           list of input BAM filenames, one per line [null]
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      <int>            read length threshold (ignore reads shorter than <int>)
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#q'
    type:
      - 'null'
      - int
    description: |
      <int>            base quality threshold
    inputBinding:
      position: 1
      prefix: '-q'
  - id: '#Q'
    type:
      - 'null'
      - int
    description: |
      <int>            mapping quality threshold
    inputBinding:
      position: 1
      prefix: '-Q'
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: |
      <chr:from-to>    region
    inputBinding:
      position: 1
      prefix: '-r'
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
  - samtools
  - depth
description: |-
  Usage: samtools depth [options] in1.bam [in2.bam [...]]
  Options:
     -b <bed>            list of positions or regions
     -f <list>           list of input BAM filenames, one per line [null]
     -l <int>            read length threshold (ignore reads shorter than <int>)
     -q <int>            base quality threshold
     -Q <int>            mapping quality threshold
     -r <chr:from-to>    region


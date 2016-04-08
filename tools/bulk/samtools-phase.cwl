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
  - id: '#inbam'
    type: File
    description: '<in.bam>'
    inputBinding:
      position: 2
  - id: '#k'
    type:
      - 'null'
      - int
    description: |
      INT    block length [13]
    inputBinding:
      position: 1
      prefix: '-k'
  - id: '#b'
    type:
      - 'null'
      - string
    description: |
      STR    prefix of BAMs to output [null]
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#q'
    type:
      - 'null'
      - int
    description: |
      INT    min het phred-LOD [37]
    inputBinding:
      position: 1
      prefix: '-q'
  - id: '#Q'
    type:
      - 'null'
      - int
    description: |
      INT    min base quality in het calling [13]
    inputBinding:
      position: 1
      prefix: '-Q'
  - id: '#D'
    type:
      - 'null'
      - int
    description: |
      INT    max read depth [256]
    inputBinding:
      position: 1
      prefix: '-D'
  - id: '#F'
    type:
      - 'null'
      - boolean
    description: "       do not attempt to fix chimeras\n"
    inputBinding:
      position: 1
      prefix: '-F'
  - id: '#A'
    type:
      - 'null'
      - boolean
    description: "       drop reads with ambiguous phase\n"
    inputBinding:
      position: 1
      prefix: '-A'
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
  - phase
description: |-
  Usage:   samtools phase [options] <in.bam>

  Options: -k INT    block length [13]
           -b STR    prefix of BAMs to output [null]
           -q INT    min het phred-LOD [37]
           -Q INT    min base quality in het calling [13]
           -D INT    max read depth [256]
           -F        do not attempt to fix chimeras
           -A        drop reads with ambiguous phase


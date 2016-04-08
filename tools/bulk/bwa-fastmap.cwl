#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: bwa-docker.cwl
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#infq'
    type: File
    description: '<in.fq>'
    inputBinding:
      position: 3
  - id: '#idxbase'
    type: boolean
    description: '<idxbase>'
    inputBinding:
      position: 2
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      INT    min SMEM length to output [17]
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#w'
    type:
      - 'null'
      - int
    description: |
      INT    max interval size to find coordiantes [20]
    inputBinding:
      position: 1
      prefix: '-w'
  - id: '#i'
    type:
      - 'null'
      - int
    description: |
      INT    min SMEM interval size [1]
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      INT    max MEM length [2147483647]
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#I'
    type:
      - 'null'
      - int
    description: |
      INT    stop if MEM is longer than -l with a size less than INT [0]
    inputBinding:
      position: 1
      prefix: '-I'
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
  - bwa
  - fastmap
description: |-
  Usage:   bwa fastmap [options] <idxbase> <in.fq>

  Options: -l INT    min SMEM length to output [17]
           -w INT    max interval size to find coordiantes [20]
           -i INT    min SMEM interval size [1]
           -l INT    max MEM length [2147483647]
           -I INT    stop if MEM is longer than -l with a size less than INT [0]


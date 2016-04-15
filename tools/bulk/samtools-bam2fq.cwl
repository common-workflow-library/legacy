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
      position: 4
  - id: '#n'
    type:
      - 'null'
      - boolean
    description: |
      don't append /1 and /2 to the read name
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#O'
    type:
      - 'null'
      - boolean
    description: "       output quality in the OQ tag if present\n"
    inputBinding:
      position: 1
      prefix: '-O'
  - id: '#s'
    type:
      - 'null'
      - File
    description: |
      FILE   write singleton reads to FILE [assume single-end]
    inputBinding:
      position: 1
      prefix: '-s'
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
  - bam2fq
description: |-
  Usage:   samtools bam2fq [-nO] [-s <outSE.fq>] <in.bam>

  Options: -n        don't append /1 and /2 to the read name
           -O        output quality in the OQ tag if present
           -s FILE   write singleton reads to FILE [assume single-end]


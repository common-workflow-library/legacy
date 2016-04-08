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
  - id: '#[out]'
    type:
      - 'null'
      - boolean
    description: '[out.]'
    inputBinding:
      position: 5
  - id: '#inbam'
    type: File
    description: '<in.bam>'
    inputBinding:
      position: 4
  - id: '#b'
    type:
      - 'null'
      - boolean
    description: "      Generate BAI-format index for BAM files [default]\n"
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: "      Generate CSI-format index for BAM files\n"
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#m'
    type:
      - 'null'
      - int
    description: |
      INT   Set minimum interval size for CSI indices to 2^INT [14]
    inputBinding:
      position: 1
      prefix: '-m'
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
  - index
description: |-
  Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
  Options:
    -b       Generate BAI-format index for BAM files [default]
    -c       Generate CSI-format index for BAM files
    -m INT   Set minimum interval size for CSI indices to 2^INT [14]


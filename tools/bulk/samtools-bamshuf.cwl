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
  - id: '#outprefix'
    type: File
    description: '<out.prefix>'
    inputBinding:
      position: 6
  - id: '#inbam'
    type: File
    description: '<in.bam>'
    inputBinding:
      position: 5
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: cLevel
    inputBinding:
      position: 4
      prefix: '-c'
  - id: '#O'
    type:
      - 'null'
      - boolean
    description: |
      output to stdout
    inputBinding:
      position: 1
      prefix: '-O'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: "     uncompressed BAM output\n"
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      INT  compression level [1]
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#n'
    type:
      - 'null'
      - int
    description: |
      INT  number of temporary files [64]
    inputBinding:
      position: 1
      prefix: '-n'
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
  - bamshuf
description: |-
  Usage:   samtools bamshuf [-Ou] [-n nFiles] [-c cLevel] <in.bam> <out.prefix>

  Options: -O      output to stdout
           -u      uncompressed BAM output
           -l INT  compression level [1]
           -n INT  number of temporary files [64]


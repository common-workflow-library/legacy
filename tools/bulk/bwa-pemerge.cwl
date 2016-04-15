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
  - id: '#[read2fq]'
    type:
      - 'null'
      - File
    description: '[read2.fq]'
    inputBinding:
      position: 4
  - id: '#read1fq'
    type: File
    description: '<read1.fq>'
    inputBinding:
      position: 3
  - id: '#m'
    type:
      - 'null'
      - boolean
    description: |
      output merged reads only
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: "      output unmerged reads only\n"
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#t'
    type:
      - 'null'
      - int
    description: |
      INT   number of threads [1]
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#T'
    type:
      - 'null'
      - int
    description: |
      INT   minimum end overlap [10]
    inputBinding:
      position: 1
      prefix: '-T'
  - id: '#Q'
    type:
      - 'null'
      - int
    description: |
      INT   max sum of errors [70]
    inputBinding:
      position: 1
      prefix: '-Q'
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
  - pemerge
description: |-
  Usage:   bwa pemerge [-mu] <read1.fq> [read2.fq]

  Options: -m       output merged reads only
           -u       output unmerged reads only
           -t INT   number of threads [1]
           -T INT   minimum end overlap [10]
           -Q INT   max sum of errors [70]


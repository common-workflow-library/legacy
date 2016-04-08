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
  - id: '#[idxbase]'
    type:
      - 'null'
      - boolean
    description: '[idxbase]'
    inputBinding:
      position: 4
  - id: '#[-d-l]'
    type:
      - 'null'
      - boolean
    description: '[-d|-l]'
    inputBinding:
      position: 2
  - id: '#d'
    type:
      - 'null'
      - boolean
    description: |
      destroy all indices in shared memory
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#l'
    type:
      - 'null'
      - boolean
    description: "      list names of indices in shared memory\n"
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#f'
    type:
      - 'null'
      - File
    description: |
      FILE  temporary file to reduce peak memory
    inputBinding:
      position: 1
      prefix: '-f'
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
  - shm
description: |-
  Usage: bwa shm [-d|-l] [-f tmpFile] [idxbase]

  Options: -d       destroy all indices in shared memory
           -l       list names of indices in shared memory
           -f FILE  temporary file to reduce peak memory


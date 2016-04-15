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
  - id: '#outbwt'
    type: File
    description: '<out.bwt>'
    inputBinding:
      position: 4
  - id: '#inpac'
    type: File
    description: '<in.pac>'
    inputBinding:
      position: 3
  - id: '#d'
    type:
      - 'null'
      - boolean
    description: ''
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
  - bwa
  - pac2bwt
description: 'Usage: bwa pac2bwt [-d] <in.pac> <out.bwt>'


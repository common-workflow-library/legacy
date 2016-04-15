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
  - id: '#bwtgen'
    type: boolean
    description: bwtgen
    inputBinding:
      position: 2
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
  - pac2bwtgen
description: 'Usage: bwtgen <in.pac> <out.bwt>'


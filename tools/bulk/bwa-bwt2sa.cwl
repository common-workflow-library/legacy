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
  - id: '#outsa'
    type: File
    description: '<out.sa>'
    inputBinding:
      position: 4
  - id: '#inbwt'
    type: File
    description: '<in.bwt>'
    inputBinding:
      position: 3
  - id: '#i'
    type:
      - 'null'
      - int
    description: '32'
    inputBinding:
      position: 2
      prefix: '-i'
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
  - bwt2sa
description: 'Usage: bwa bwt2sa [-i 32] <in.bwt> <out.sa>'


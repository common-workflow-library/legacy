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
  - id: '#[outprefix]'
    type:
      - 'null'
      - File
    description: '[<out.prefix>]'
    inputBinding:
      position: 4
  - id: '#infasta'
    type: File
    description: '<in.fasta>'
    inputBinding:
      position: 3
  - id: '#f'
    type:
      - 'null'
      - boolean
    description: ''
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
  - fa2pac
description: 'Usage: bwa fa2pac [-f] <in.fasta> [<out.prefix>]'


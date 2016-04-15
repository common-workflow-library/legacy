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
      position: 7
  - id: '#insai'
    type: File
    description: '<in.sai>'
    inputBinding:
      position: 6
  - id: '#prefix'
    type: File
    description: '<prefix>'
    inputBinding:
      position: 5
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: RG_line
    inputBinding:
      position: 4
      prefix: '-r'
  - id: '#f'
    type:
      - 'null'
      - File
    description: out.sam
    inputBinding:
      position: 3
      prefix: '-f'
  - id: '#n'
    type:
      - 'null'
      - boolean
    description: max_occ
    inputBinding:
      position: 2
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
  - bwa
  - samse
description: 'Usage: bwa samse [-n max_occ] [-f out.sam] [-r RG_line] <prefix> <in.sai> <in.fq>'


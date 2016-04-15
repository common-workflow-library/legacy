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
      position: 8
  - id: '#f'
    type:
      - 'null'
      - boolean
    description: ref
    inputBinding:
      position: 7
      prefix: '-f'
  - id: '#2'
    type:
      - 'null'
      - boolean
    description: em2
    inputBinding:
      position: 6
      prefix: '-2'
  - id: '#1'
    type:
      - 'null'
      - boolean
    description: em1
    inputBinding:
      position: 5
      prefix: '-1'
  - id: '#0'
    type:
      - 'null'
      - boolean
    description: em0
    inputBinding:
      position: 4
      prefix: '-0'
  - id: '#i'
    type:
      - 'null'
      - boolean
    description: inPen
    inputBinding:
      position: 3
      prefix: '-i'
  - id: '#Q'
    type:
      - 'null'
      - boolean
    description: minQ
    inputBinding:
      position: 2
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
  - samtools
  - targetcut
description: 'Usage: samtools targetcut [-Q minQ] [-i inPen] [-0 em0] [-1 em1] [-2 em2] [-f ref] <in.bam>'


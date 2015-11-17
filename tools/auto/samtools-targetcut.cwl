#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: samtools-docker.cwl
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#inbam'
    type: File
    description: '<in.bam>'
    inputBinding:
      position: 8
  - id: '#f'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#2'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-2'
  - id: '#1'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-1'
  - id: '#0'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-0'
  - id: '#i'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#Q'
    type: boolean
    description: ''
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
  - samtools
  - targetcut
description: 'Usage: samtools targetcut [-Q minQ] [-i inPen] [-0 em0] [-1 em1] [-2 em2] [-f ref] <in.bam>'


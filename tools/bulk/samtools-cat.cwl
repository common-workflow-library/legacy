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
  - id: '#[]'
    type:
      - 'null'
      - boolean
    description: '[...]'
    inputBinding:
      position: 6
  - id: '#in2bam'
    type: File
    description: '<in2.bam>'
    inputBinding:
      position: 5
  - id: '#in1bam'
    type: File
    description: '<in1.bam>'
    inputBinding:
      position: 4
  - id: '#o'
    type:
      - 'null'
      - File
    description: out.bam
    inputBinding:
      position: 3
      prefix: '-o'
  - id: '#h'
    type:
      - 'null'
      - File
    description: header.sam
    inputBinding:
      position: 2
      prefix: '-h'
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
  - cat
description: 'Usage: samtools cat [-h header.sam] [-o out.bam] <in1.bam> <in2.bam> [...]'


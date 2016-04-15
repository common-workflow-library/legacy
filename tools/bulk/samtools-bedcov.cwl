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
      position: 4
  - id: '#in1bam'
    type: File
    description: '<in1.bam>'
    inputBinding:
      position: 3
  - id: '#inbed'
    type: File
    description: '<in.bed>'
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
  - samtools
  - bedcov
description: 'Usage: samtools bedcov <in.bed> <in1.bam> [...]'


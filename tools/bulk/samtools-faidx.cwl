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
  - id: '#]'
    type: boolean
    description: ']'
    inputBinding:
      position: 4
  - id: '#filefafilefagz'
    type: File
    description: '<file.fa|file.fa.gz>'
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
  - faidx
description: 'Usage:   samtools faidx <file.fa|file.fa.gz> [<reg> [...]]'


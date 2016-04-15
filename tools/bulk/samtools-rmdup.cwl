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
  - id: '#outputbam'
    type: File
    description: '<output.bam>'
    inputBinding:
      position: 4
  - id: '#inputsrtbam'
    type: File
    description: '<input.srt.bam>'
    inputBinding:
      position: 3
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      rmdup for SE reads
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: "   treat PE reads as SE in rmdup (force -s)\n"
    inputBinding:
      position: 1
      prefix: '-S'
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
  - rmdup
description: |-
  Usage:  samtools rmdup [-sS] <input.srt.bam> <output.bam>

  Option: -s    rmdup for SE reads
          -S    treat PE reads as SE in rmdup (force -s)


#!/usr/bin/env cwl-runner
#
# To use it as stand alone tool current working directory should not have input .fa file
#    example: "./samtools-faidx.cwl --input=./test-files/mm10.fa"
# As part of a workflow should be no problem at all
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: samtools-docker.cwl
  - class: CreateFileRequirement
    fileDef:
      - filename:
          engine: node-engine.cwl
          script: $job['input'].path.split('/').slice(-1)[0]
        fileContent:
          engine: "cwl:JsonPointer"
          script: /job/input
inputs:
  - id: '#input'
    type: File
    description: '<file.fa|file.fa.gz>'

  - id: '#filename'
    type: ["null",string]
    default: ""
    inputBinding:
      position: 1
      valueFrom:
        engine: node-engine.cwl
        script: $job['input'].path.split('/').slice(-1)[0]

  - id: '#region'
    type: ["null",string]
    inputBinding:
      position: 2

outputs:
  - id: "#index_result"
    type: File
    outputBinding:
      glob:
        engine: node-engine.cwl
        script: $job['input'].path.split('/').slice(-1)[0]+'.fai'

baseCommand:
  - samtools
  - faidx
description: 'Usage:   samtools faidx <file.fa|file.fa.gz> [<reg> [...]]'


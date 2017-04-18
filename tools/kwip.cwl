#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

hints:
  - class: ResourceRequirement
    coresMin: 4

inputs:
  - id: unweighted
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --unweighted
  - id: hashes
    type:
      type: array
      items: File
    inputBinding:
      position: 2
  - id: distance
    type: string
    inputBinding:
      prefix: --distance
  - id: kernel
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --kernel

outputs:
  - id: kernelmat
    type:
      - 'null'
      - File
    outputBinding:
      glob: '*.kern'
  - id: distancemat
    type: File
    outputBinding:
      glob: '*.dist'

baseCommand: [kwip,]

arguments:
  - valueFrom: $(runtime.cores)
    position: 1
    prefix: --threads

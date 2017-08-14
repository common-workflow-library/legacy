#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

inputs:
  files:
    type: File[]
    inputBinding: {position: 1}
    streamable: true

baseCommand: [wc, -l]

outputs:
  counts: stdout

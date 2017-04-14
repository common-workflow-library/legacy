#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

inputs:
  infiles:
    type: File[]
    inputBinding: {position: 1}

outputs:
  outfile:
    type: stdout

baseCommand: [wc, -l]


#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  aligned:
    type: File
    inputBinding:
      prefix: "-sequence"

  method:
    type:
      type: enum
      symbols: [ f, j, k ]
    inputBinding:
      prefix: "-method"

baseCommand: [ fnadist ]

arguments: [ "-auto", aligned.dend, "-outfile", aligned.fdnadist ]

outputs:
  distance:
    type: File
    outputBinding:
      glob: aligned.fdnadist

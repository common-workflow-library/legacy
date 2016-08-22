#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  dist:
    type: File
    inputBinding:
      prefix: "-datafile"

baseCommand: [ fneighbor ]

arguments: [ "-outfile", aligned.tree ]

outputs:
  result:
    type: File
    outputBinding:
      glob: aligned.tree

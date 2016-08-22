#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  unaligned:
    type: File
    inputBinding:
      prefix: "-sequence"

baseCommand: [ emma ]

arguments: [ "-dendoutfile", aligned.dend, "-outseq", "aligned.fa" ]

outputs:
  alignedseq:
    type: File
    outputBinding:
      glob: aligned.fa

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  fastqFile:
    type: File
    inputBinding:
      position: 1

baseCommand: [ fastqc, "--outdir", .] #"--extract"]

outputs:
  zippedFile:
    type: File
    outputBinding:
      glob: "*.zip"

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- $import: fastqc_dockerfile.yml

inputs:
  fastqFile:
    type: File #No reason to aaccept multiple files as no overall report is generated
    inputBinding:
      position: 1

baseCommand: [ fastqc, "--outdir", . , "--extract"]

outputs:
  zippedFile:
    type: File
    outputBinding:
      glob: "*.zip"
  report:
    type: Directory
    outputBinding:
      glob: "."

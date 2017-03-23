cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement

inputs:
  infile:
    type: File
    inputBinding:
      position: 1

outputs:
  outfile:
    type: File
    outputSource: complement/outfile

steps:
  reverse:
    run: reverse.cwl
    in:
      infile: infile
    out: [outfile]
  complement:
    run: complement.cwl
    in:
      infile: reverse/outfile
    out: [outfile]

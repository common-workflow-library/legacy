#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement

inputs:
  infile:
    type: File
    inputBinding:
      position: 1
  outfile_name:
    type: string?

outputs:
  revcomp_dnafile:
    type: File
    outputSource: rename/outfile

steps:
  reverse:
    run: reverse.cwl
    in:
      dnafile: infile
    out: [rev_dnafile]
  complement:
    run: complement.cwl
    in:
      dnafile: reverse/rev_dnafile
    out: [comp_dnafile]
  rename:
    run:
      class: ExpressionTool
      inputs:
        infile:
          type: File
        outfile_name:
          type: string?
      outputs:
        outfile: File
      expression: >
        ${
        var outfile = inputs.infile;
        if (inputs.outfile_name) {
          outfile.basename = inputs.outfile_name;
        }
        return { "outfile": outfile }; }
    in:
      infile: complement/comp_dnafile
      outfile_name: outfile_name
    out: [outfile]

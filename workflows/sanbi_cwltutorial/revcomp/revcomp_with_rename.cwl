cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement

inputs:
  infile:
    type: File
    inputBinding:
      position: 1
  outfile:
    type: string?

outputs:
  outfile:
    type: File
    outputSource: rename/outfile

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
  rename:
    run:
      class: ExpressionTool
      inputs:
        infile:
          type: File
        outfile:
          type: string?
      outputs:
        outfile: File
      expression: >
        ${
        var outfile = inputs.infile;
        if (inputs.outfile) {
          outfile.basename = inputs.outfile;
        }
        return { outfile: outfile }; }
    in:
      infile: complement/outfile
      outfile: outfile
    out: [outfile]

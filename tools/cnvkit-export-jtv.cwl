#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'export', 'jtv']

description: |
  Convert log2 ratios to Java TreeView's native format.

inputs:
  
  filenames:
    type:
      type: array
      items: string
  
    description: Log2 copy ratio data file(s) (*.cnr), the output of the 'fix' sub-command.
    inputBinding:
      position: 1

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

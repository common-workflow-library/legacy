#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'gender']

description: |
  Guess samples' gender from the relative coverage of chromosome X.

inputs:
  
  targets:
    type:
      type: array
      items: string
  
    description: Copy number or copy ratio files (*.cnn, *.cnr).
    inputBinding:
      position: 1

  male_reference:
    type: ["null", boolean]
    default: False
    description: Assume inputs are already normalized to a male reference (i.e. female samples will have +1 log-coverage of chrX; otherwise male samples would have -1 chrX).
    inputBinding:
      prefix: --male-reference 

  output:
    type: ["null", string]
    description: Output table file name.
    inputBinding:
      prefix: --output 


outputs:
    []

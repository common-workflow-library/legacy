#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'export', 'seg']

description: |
  Convert segments to SEG format.
  
      Compatible with IGV and GenePattern.
      

inputs:
  
  filenames:
    type:
      type: array
      items: string
  
    description: Segmented copy ratio data file(s) (*.cns), the output of the 'segment' sub-command.
    inputBinding:
      position: 1

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

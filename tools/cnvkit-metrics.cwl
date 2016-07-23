#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'metrics']

description: |
  Compute coverage deviations and other metrics for self-evaluation.
      

inputs:
  
  cnarrays:
    type:
      type: array
      items: string
  
    description: One or more bin-level coverage data files (*.cnn, *.cnr).
    inputBinding:
      position: 1

  segments:
    type:
      type: array
      items: string
  
    description: One or more segmentation data files (*.cns, output of the 'segment' command). If more than one file is given, the number must match the coverage data files, in which case the input files will be paired together in the given order. Otherwise, the same segments will be used for all coverage files.
    inputBinding:
      prefix: --segments 

  output:
    type: ["null", string]
    description: Output table file name.
    inputBinding:
      prefix: --output 


outputs:
    []

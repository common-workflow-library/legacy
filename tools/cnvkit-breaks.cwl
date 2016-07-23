#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'breaks']

description: |
  List the targeted genes in which a copy number breakpoint occurs.

inputs:
  
  filename:
    type: string
  
    description: Processed sample coverage data file (*.cnr), the output of the 'fix' sub-command.
    inputBinding:
      position: 1

  segment:
    type: string
  
    description: Segmentation calls (.cns), the output of the 'segment' command).
    inputBinding:
      position: 2

  min_probes:
    type: ["null", int]
    default: 1
    description: Minimum number of within-gene probes on both sides of a breakpoint to report it. [Default - %(default)d]
    inputBinding:
      prefix: --min-probes 

  output:
    type: ["null", string]
    description: Output table file name.
    inputBinding:
      prefix: --output 


outputs:
    []

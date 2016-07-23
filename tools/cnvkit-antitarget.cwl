#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'antitarget']

description: |
  Derive a background/antitarget BED file from a target BED file.

inputs:
  
  interval:
    type: string
  
    description: BED or interval file listing the targeted regions.
    inputBinding:
      position: 1

  access:
    type: ["null", string]
    description: Regions of accessible sequence on chromosomes (.bed), as output by genome2access.py.
    inputBinding:
      prefix: --access 

  avg_size:
    type: ["null", int]
    default: 100000
    description: Average size of antitarget bins (results are approximate). [Default - %(default)s]
    inputBinding:
      prefix: --avg-size 

  min_size:
    type: ["null", int]
    description: Minimum size of antitarget bins (smaller regions are dropped). [Default - 1/16 avg size, calculated]
    inputBinding:
      prefix: --min-size 

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

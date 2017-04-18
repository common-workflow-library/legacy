#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'target']

description: |
  Transform bait intervals into targets more suitable for CNVkit.

inputs:
  
  interval:
    type: string
  
    description: BED or interval file listing the targeted regions.
    inputBinding:
      position: 1

  annotate:
    type: ["null", string]
    description: UCSC refFlat.txt or ensFlat.txt file for the reference genome. Pull gene names from this file and assign them to the target regions.
    inputBinding:
      prefix: --annotate 

  short_names:
    type: ["null", boolean]
    default: False
    description: Reduce multi-accession bait labels to be short and consistent.
    inputBinding:
      prefix: --short-names 

  split:
    type: ["null", boolean]
    default: False
    description: Split large tiled intervals into smaller, consecutive targets.
    inputBinding:
      prefix: --split 

  avg_size:
    type: ["null", int]
    default: 266.6666666666667
    description: Average size of split target bins (results are approximate). [Default - %(default)s]
    inputBinding:
      prefix: --avg-size 

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

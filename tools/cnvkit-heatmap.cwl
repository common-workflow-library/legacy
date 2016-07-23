#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'heatmap']

description: |
  Plot copy number for multiple samples as a heatmap.

inputs:
  
  filenames:
    type:
      type: array
      items: string
  
    description: Sample coverages as raw probes (.cnr) or segments (.cns).
    inputBinding:
      position: 1

  chromosome:
    type: ["null", string]
    description: Chromosome (e.g. 'chr1') or chromosomal range (e.g. 'chr1 -2333000-2444000') to display. If a range is given, all targeted genes in this range will be shown, unless '--gene'/'-g' is already given.
    inputBinding:
      prefix: --chromosome 

  desaturate:
    type: ["null", boolean]
    default: False
    description: Tweak color saturation to focus on significant changes.
    inputBinding:
      prefix: --desaturate 

  output:
    type: ["null", string]
    description: Output PDF file name.
    inputBinding:
      prefix: --output 


outputs:
    []

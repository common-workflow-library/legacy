#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'diagram']

description: |
  Draw copy number (log2 coverages, CBS calls) on chromosomes as a diagram.
  
      If both the raw probes and segments are given, show them side-by-side on
      each chromosome (segments on the left side, probes on the right side).
      

inputs:
  
  filename:
    type: ["null", string]
    description: Processed coverage data file (*.cnr), the output of the 'fix' sub-command.
    inputBinding:
      position: 1

  segment:
    type: ["null", string]
    description: Segmentation calls (.cns), the output of the 'segment' command.
    inputBinding:
      prefix: --segment 

  threshold:
    type: ["null", float]
    default: 0.5
    description: Copy number change threshold to label genes. [Default - %(default)s]
    inputBinding:
      prefix: --threshold 

  min_probes:
    type: ["null", int]
    default: 3
    description: Minimum number of covered probes to label a gene. [Default - %(default)d]
    inputBinding:
      prefix: --min-probes 

  male_reference:
    type: ["null", boolean]
    default: False
    description: Assume inputs are already corrected against a male reference (i.e. female samples will have +1 log-CNR of chrX; otherwise male samples would have -1 chrX).
    inputBinding:
      prefix: --male-reference 

  output:
    type: ["null", string]
    description: Output PDF file name.
    inputBinding:
      prefix: --output 


outputs:
    []

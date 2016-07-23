#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'gainloss']

description: |
  Identify targeted genes with copy number gain or loss.

inputs:
  
  filename:
    type: string
  
    description: Processed sample coverage data file (*.cnr), the output of the 'fix' sub-command.
    inputBinding:
      position: 1

  segment:
    type: ["null", string]
    description: Segmentation calls (.cns), the output of the 'segment' command).
    inputBinding:
      prefix: --segment 

  threshold:
    type: ["null", float]
    default: 0.2
    description: Copy number change threshold to report a gene gain/loss. [Default - %(default)s]
    inputBinding:
      prefix: --threshold 

  min_probes:
    type: ["null", int]
    default: 3
    description: Minimum number of covered probes to report a gain/loss. [Default - %(default)d]
    inputBinding:
      prefix: --min-probes 

  drop_low_coverage:
    type: ["null", boolean]
    default: False
    description: Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples.
    inputBinding:
      prefix: --drop-low-coverage 

  male_reference:
    type: ["null", boolean]
    default: False
    description: Assume inputs are already corrected against a male reference (i.e. female samples will have +1 log-coverage of chrX; otherwise male samples would have -1 chrX).
    inputBinding:
      prefix: --male-reference 

  output:
    type: ["null", string]
    description: Output table file name.
    inputBinding:
      prefix: --output 


outputs:
    []

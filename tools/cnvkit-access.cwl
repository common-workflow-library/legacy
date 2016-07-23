#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'access']

description: |
  List the locations of accessible sequence regions in a FASTA file.

inputs:
  
  fa_fname:
    type: string
  
    description: Genome FASTA file name
    inputBinding:
      position: 1

  min_gap_size:
    type: ["null", int]
    default: 5000
    description: Minimum gap size between accessible sequence regions. Regions separated by less than this distance will be joined together. [Default - %(default)s]
    inputBinding:
      prefix: --min-gap-size 

  exclude:
    type:
    - "null"
    - type: array
      items: string
  
    default: []
    description: Additional regions to exclude, in BED format. Can be used multiple times.
    inputBinding:
      prefix: --exclude 

  output:
    type: ["null", File]
    description: Output file name
    inputBinding:
      prefix: --output 


outputs:

  output_out:
    type: File
  
    description: Output file name
    outputBinding:
      glob: $(inputs.output.path)

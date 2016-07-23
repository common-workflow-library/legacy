#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'import-theta']

description: |
  Convert THetA output to a BED-like, CNVkit-like tabular format.
  
      Equivalently, use the THetA results file to convert CNVkit .cns segments to
      integer copy number calls.
      

inputs:
  
  tumor_cns:
    type: string
  
  
    inputBinding:
      position: 1

  theta_results:
    type: string
  
  
    inputBinding:
      position: 2

  ploidy:
    type: ["null", int]
    default: 2
    description: Ploidy of normal cells. [Default - %(default)d]
    inputBinding:
      prefix: --ploidy 

  output_dir:
    type: ["null", string]
    default: .
    description: Output directory name.
    inputBinding:
      prefix: --output-dir 


outputs:
    []

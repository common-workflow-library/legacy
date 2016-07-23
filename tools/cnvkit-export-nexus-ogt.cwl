#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'export', 'nexus-ogt']

description: |
  Convert log2 ratios and b-allele freqs to Nexus "Custom-OGT" format.

inputs:
  
  filename:
    type: string
  
    description: Log2 copy ratio data file (*.cnr), the output of the 'fix' sub-command.
    inputBinding:
      position: 1

  vcf:
    type: string
  
    description: VCF of SNVs for the same sample, to calculate b-allele frequencies.
    inputBinding:
      position: 2

  sample_id:
    type: ["null", string]
    description: Specify the name of the sample in the VCF to use to extract b-allele frequencies.
    inputBinding:
      prefix: --sample-id 

  min_variant_depth:
    type: ["null", int]
    default: 20
    description: Minimum read depth for a SNV to be included in the b-allele frequency calculation. [Default - %(default)s]
    inputBinding:
      prefix: --min-variant-depth 

  min_weight:
    type: ["null", float]
    default: 0.0
    description: Minimum weight (between 0 and 1) for a bin to be included in the output. [Default - %(default)s]
    inputBinding:
      prefix: --min-weight 

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

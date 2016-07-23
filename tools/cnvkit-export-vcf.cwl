#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'export', 'vcf']

description: |
  Convert segments to VCF format.
  
      Input is a segmentation file (.cns) where, preferably, log2 ratios have
      already been adjusted to integer absolute values using the 'call' command.
      

inputs:
  
  segments:
    type: string
  
    description: Segmented copy ratio data file (*.cns), the output of the 'segment' or 'call' sub-commands.
    inputBinding:
      position: 1

  sample_id:
    type: ["null", string]
    description: Sample name to write in the genotype field of the output VCF file. [Default - use the sample ID, taken from the file name]
    inputBinding:
      prefix: --sample-id 

  ploidy:
    type: ["null", int]
    default: 2
    description: Ploidy of the sample cells. [Default - %(default)d]
    inputBinding:
      prefix: --ploidy 

  gender:
    type:
    - "null"
    - type: enum
      symbols: ['m', 'male', 'Male', 'f', 'female', 'Female']
    description: Specify the sample's gender as male or female. (Otherwise guessed from chrX copy number).
    inputBinding:
      prefix: --gender 

  male_reference:
    type: ["null", boolean]
    default: False
    description: Was a male reference used? If so, expect half ploidy on chrX and chrY; otherwise, only chrY has half ploidy. In CNVkit, if a male reference was used, the "neutral" copy number (ploidy) of chrX is 1; chrY is haploid for either gender reference.
    inputBinding:
      prefix: --male-reference 

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

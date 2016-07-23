#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'rescale']

description: |
  [DEPRECATED] Rescale segment copy ratios given known purity and ploidy.
  
      Instead, use the command "call -m none".
      

inputs:
  
  filename:
    type: string
  
    description: Copy ratios (.cnr or .cns).
    inputBinding:
      position: 1

  center:
    type:
    - "null"
    - type: enum
      symbols: ['mean', 'median', 'mode', 'biweight']
    description: Re-center the log2 ratio values using this estimate of the center or average value.
    inputBinding:
      prefix: --center 

  ploidy:
    type: ["null", int]
    default: 2
    description: Ploidy of the sample cells. [Default - %(default)d]
    inputBinding:
      prefix: --ploidy 

  purity:
    type: ["null", float]
    description: Estimated tumor cell fraction, a.k.a. purity or cellularity.
    inputBinding:
      prefix: --purity 

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
    description: Output table file name (CNR-like table of segments, .cns).
    inputBinding:
      prefix: --output 


outputs:
    []

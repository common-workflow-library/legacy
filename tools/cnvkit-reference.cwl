#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'reference']

description: |
  Compile a coverage reference from the given files (normal samples).

inputs:
  
  references:
    type:
    - "null"
    - type: array
      items: string
  
    description: Normal-sample target or antitarget .cnn files, or the directory that contains them.
    inputBinding:
      position: 1

  fasta:
    type: ["null", string]
    description: Reference genome, FASTA format (e.g. UCSC hg19.fa)
    inputBinding:
      prefix: --fasta 

  targets:
    type: ["null", string]
    description: Target intervals (.bed or .list)
    inputBinding:
      prefix: --targets 

  antitargets:
    type: ["null", string]
    description: Antitarget intervals (.bed or .list)
    inputBinding:
      prefix: --antitargets 

  male_reference:
    type: ["null", boolean]
    default: False
    description: Create a male reference - shift female samples' chrX log-coverage by -1, so the reference chrX average is -1. Otherwise, shift male samples' chrX by +1, so the reference chrX average is 0.
    inputBinding:
      prefix: --male-reference 

  do_gc:
    type: ["null", boolean]
    default: True
    description: Skip GC correction.
    inputBinding:
      prefix: --no-gc 

  do_edge:
    type: ["null", boolean]
    default: True
    description: Skip edge-effect correction.
    inputBinding:
      prefix: --no-edge 

  do_rmask:
    type: ["null", boolean]
    default: True
    description: Skip RepeatMasker correction.
    inputBinding:
      prefix: --no-rmask 

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

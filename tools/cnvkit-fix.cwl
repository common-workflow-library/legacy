#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'fix']

description: |
  Combine target and antitarget coverages and correct for biases.
  
      Adjust raw coverage data according to the given reference, correct potential
      biases and re-center.
      

inputs:
  
  target:
    type: string
  
    description: Target coverage file (.targetcoverage.cnn).
    inputBinding:
      position: 1

  antitarget:
    type: string
  
    description: Antitarget coverage file (.antitargetcoverage.cnn).
    inputBinding:
      position: 2

  reference:
    type: string
  
    description: Reference coverage (.cnn).
    inputBinding:
      position: 3

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

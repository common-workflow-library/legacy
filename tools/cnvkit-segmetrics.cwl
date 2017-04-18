#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'segmetrics']

description: |
  Compute segment-level metrics from bin-level log2 ratios.

inputs:
  
  cnarray:
    type: string
  
    description: Bin-level copy ratio data file (*.cnn, *.cnr).
    inputBinding:
      position: 1

  segments:
    type: string
  
    description: Segmentation data file (*.cns, output of the 'segment' command).
    inputBinding:
      prefix: --segments 

  drop_low_coverage:
    type: ["null", boolean]
    default: False
    description: Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples.
    inputBinding:
      prefix: --drop-low-coverage 

  output:
    type: ["null", string]
    description: Output table file name.
    inputBinding:
      prefix: --output 

  stdev:
    type: ["null", boolean]
    default: False
    description: Standard deviation.
    inputBinding:
      prefix: --stdev 

  mad:
    type: ["null", boolean]
    default: False
    description: Median absolute deviation (standardized).
    inputBinding:
      prefix: --mad 

  iqr:
    type: ["null", boolean]
    default: False
    description: Inter-quartile range.
    inputBinding:
      prefix: --iqr 

  bivar:
    type: ["null", boolean]
    default: False
    description: Tukey's biweight midvariance.
    inputBinding:
      prefix: --bivar 

  ci:
    type: ["null", boolean]
    default: False
    description: Confidence interval (by bootstrap).
    inputBinding:
      prefix: --ci 

  pi:
    type: ["null", boolean]
    default: False
    description: Prediction interval.
    inputBinding:
      prefix: --pi 

  alpha:
    type: ["null", float]
    default: 0.05
    description: Level to estimate confidence and prediction intervals; use with --ci and --pi. [Default - %(default)s]
    inputBinding:
      prefix: --alpha 

  bootstrap:
    type: ["null", int]
    default: 100
    description: Number of bootstrap iterations to estimate confidence interval; use with --ci. [Default - %(default)d]
    inputBinding:
      prefix: --bootstrap 


outputs:
    []

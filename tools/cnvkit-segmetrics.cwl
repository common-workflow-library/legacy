#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'segmetrics']

requirements:
  - class: InlineJavascriptRequirement

description: |
  Compute segment-level metrics from bin-level log2 ratios.

inputs:
  

- id: cnarray
  type: string

  description: Bin-level copy ratio data file (*.cnn, *.cnr).
  inputBinding:
    position: 1

- id: segments
  type: string

  description: Segmentation data file (*.cns, output of the 'segment' command).
  inputBinding:
    prefix: --segments 

- id: drop_low_coverage
  type: ["null", boolean]
  default: null
  description: Drop very-low-coverage bins before segmentation to avoid
                false-positive deletions in poor-quality tumor samples.
  inputBinding:
    prefix: --drop-low-coverage 

- id: output
  type: ["null", string]
  description: Output table file name.
  inputBinding:
    prefix: --output 

- id: stdev
  type: ["null", boolean]
  default: null
  description: Standard deviation.
  inputBinding:
    prefix: --stdev 

- id: mad
  type: ["null", boolean]
  default: null
  description: Median absolute deviation (standardized).
  inputBinding:
    prefix: --mad 

- id: iqr
  type: ["null", boolean]
  default: null
  description: Inter-quartile range.
  inputBinding:
    prefix: --iqr 

- id: bivar
  type: ["null", boolean]
  default: null
  description: Tukey's biweight midvariance.
  inputBinding:
    prefix: --bivar 

- id: ci
  type: ["null", boolean]
  default: null
  description: Confidence interval (by bootstrap).
  inputBinding:
    prefix: --ci 

- id: pi
  type: ["null", boolean]
  default: null
  description: Prediction interval.
  inputBinding:
    prefix: --pi 

- id: alpha
  type: ["null", float]
  default: 0.05
  description: Level to estimate confidence and prediction intervals;
                use with --ci and --pi. [Default - %(default)s]
  inputBinding:
    prefix: --alpha 

- id: bootstrap
  type: ["null", int]
  default: 100
  description: Number of bootstrap iterations to estimate confidence interval;
                use with --ci. [Default - %(default)d]
  inputBinding:
    prefix: --bootstrap 

outputs:
    []

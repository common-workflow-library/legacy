#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'segment']

description: |
  Infer copy number segments from the given coverage table.

inputs:
  
  filename:
    type: string
  
    description: Bin-level log2 ratios (.cnr file), as produced by 'fix'.
    inputBinding:
      position: 1

  output:
    type: ["null", string]
    description: Output table file name (CNR-like table of segments, .cns).
    inputBinding:
      prefix: --output 

  dataframe:
    type: ["null", string]
    description: File name to save the raw R dataframe emitted by CBS or Fused Lasso. (Useful for debugging.)
    inputBinding:
      prefix: --dataframe 

  method:
    type:
    - "null"
    - type: enum
      symbols: ['cbs', 'haar', 'flasso']
    default: cbs
    description: Segmentation method (CBS, HaarSeg, or Fused Lasso). [Default - %(default)s]
    inputBinding:
      prefix: --method 

  threshold:
    type: ["null", float]
    description: Significance threshold (p-value or FDR, depending on method) to accept breakpoints during segmentation.
    inputBinding:
      prefix: --threshold 

  vcf:
    type: ["null", string]
    description: VCF file name containing variants for segmentation by allele frequencies.
    inputBinding:
      prefix: --vcf 

  drop_low_coverage:
    type: ["null", boolean]
    default: False
    description: Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples.
    inputBinding:
      prefix: --drop-low-coverage 

  drop_outliers:
    type: ["null", float]
    default: 10
    description: Drop outlier bins more than this many multiples of the 95th quantile away from the average within a rolling window. Set to 0 for no outlier filtering. [Default - %(default)g]
    inputBinding:
      prefix: --drop-outliers 

  rlibpath:
    type: ["null", string]
    description: Path to an alternative site-library to use for R packages.
    inputBinding:
      prefix: --rlibpath 


outputs:
    []

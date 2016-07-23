#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'loh']

description: |
  [DEPRECATED] Plot allelic frequencies at each variant position in a VCF file.
  
      Divergence from 0.5 indicates loss of heterozygosity in a tumor sample.
  
      Instead, use the command "scatter -v".
      

inputs:
  
  variants:
    type: string
  
    description: Sample variants in VCF format.
    inputBinding:
      position: 1

  segment:
    type: ["null", string]
    description: Segmentation calls (.cns), the output of the 'segment' command.
    inputBinding:
      prefix: --segment 

  min_depth:
    type: ["null", int]
    default: 20
    description: Minimum read depth for a variant to be displayed. [Default - %(default)s]
    inputBinding:
      prefix: --min-depth 

  sample_id:
    type: ["null", string]
    description: Sample name to use for LOH calculations from the input VCF.
    inputBinding:
      prefix: --sample-id 

  normal_id:
    type: ["null", string]
    description: Corresponding normal sample ID in the input VCF.
    inputBinding:
      prefix: --normal-id 

  trend:
    type: ["null", boolean]
    default: False
    description: Draw a smoothed local trendline on the scatter plot.
    inputBinding:
      prefix: --trend 

  output:
    type: ["null", string]
    description: Output PDF file name.
    inputBinding:
      prefix: --output 


outputs:
    []

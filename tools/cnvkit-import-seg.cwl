#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'import-seg']

description: |
  Convert a SEG file to CNVkit .cns files.

inputs:
  
  segfile:
    type: string
  
    description: Input file in SEG format. May contain multiple samples.
    inputBinding:
      position: 1

  chromosomes:
    type: ["null", string]
    description: Mapping of chromosome indexes to names. Syntax - "from1 -to1,from2 -to2". Or use "human" for the preset - "23 -X,24 -Y,25 -M".
    inputBinding:
      prefix: --chromosomes 

  prefix:
    type: ["null", string]
    description: Prefix to add to chromosome names (e.g 'chr' to rename '8' in the SEG file to 'chr8' in the output).
    inputBinding:
      prefix: --prefix 

  from_log10:
    type: ["null", boolean]
    default: False
    description: Convert base-10 logarithm values in the input to base-2 logs.
    inputBinding:
      prefix: --from-log10 

  output_dir:
    type: ["null", string]
    default: .
    description: Output directory name.
    inputBinding:
      prefix: --output-dir 


outputs:
    []

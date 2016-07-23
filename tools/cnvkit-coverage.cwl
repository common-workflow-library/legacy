#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'coverage']

description: |
  Calculate coverage in the given regions from BAM read depths.

inputs:
  
  bam_file:
    type: string
  
    description: Mapped sequence reads (.bam)
    inputBinding:
      position: 1

  interval:
    type: string
  
    description: Intervals (.bed or .list)
    inputBinding:
      position: 2

  count:
    type: ["null", boolean]
    default: False
    description: Get read depths by counting read midpoints within each bin. (An alternative algorithm).
    inputBinding:
      prefix: --count 

  min_mapq:
    type: ["null", int]
    default: 0
    description: Minimum mapping quality score (phred scale 0-60) to count a read for coverage depth. [Default - %(default)s]
    inputBinding:
      prefix: --min-mapq 

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

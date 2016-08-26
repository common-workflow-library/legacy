#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- $import: samtools-docker.yml
- class: InlineJavascriptRequirement

inputs:
  compression_level:
    type: int?
    inputBinding:
      prefix: -l
    doc: |
      Set compression level, from 0 (uncompressed) to 9 (best)
  threads:
    type: int?
    inputBinding:
      prefix: -@

    doc: Set number of sorting and compression threads [1]
  memory:
    type: string?
    inputBinding:
      prefix: -m
    doc: |
      Set maximum memory per thread; suffix K/M/G recognized [768M]
  input:
    type: File
    inputBinding:
      position: 1

    doc: Input bam file.
  #output_name:
  #  type: string
  #  inputBinding:
  #    position: 2

  output_name:
    type: string
    inputBinding:
      prefix: -o
    doc: Desired output filename.

  temp_dir:
    type: string?
    default: $TMPDIR
    inputBinding:
      position: 1
      prefix: -T
    doc: |
      working temp directory

  sort_by_name:
    type: boolean?
    inputBinding:
      prefix: -n

    doc: Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
outputs:
  sorted:
    type: File
    format: http://edamontology.org/format_2572
    outputBinding:
      glob: $(inputs.output_name)

baseCommand: [samtools, sort]

#doc: |
#  samtools-sort.cwl is developed for CWL consortium
#    Usage: samtools sort [options...] [in.bam]
#    Options:
#      -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
#      -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
#      -n         Sort by read name
#      -o FILE    Write final output to FILE rather than standard output
#      -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')   (either -O or
#      -T PREFIX  Write temporary files to PREFIX.nnnn.bam       -T is required)
#      -@ INT     Set number of sorting and compression threads [1]
#
#    Legacy usage: samtools sort [options...] <in.bam> <out.prefix>
#    Options:
#      -f         Use <out.prefix> as full final filename rather than prefix
#      -o         Write final output to stdout rather than <out.prefix>.bam
#      -l,m,n,@   Similar to corresponding options above
#
#

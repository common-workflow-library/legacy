#!/usr/bin/env cwl-runner

class: CommandLineTool

description: |
  Usage: samtools sort [options...] [in.bam]
  Options:
    -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
    -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
    -n         Sort by read name
    -o FILE    Write final output to FILE rather than standard output
    -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')   (either -O or
    -T PREFIX  Write temporary files to PREFIX.nnnn.bam       -T is required)
    -@ INT     Set number of sorting and compression threads [1]

  Legacy usage: samtools sort [options...] <in.bam> <out.prefix>
  Options:
    -f         Use <out.prefix> as full final filename rather than prefix
    -o         Write final output to stdout rather than <out.prefix>.bam
    -l,m,n,@   Similar to corresponding options above

requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: samtools-docker.cwl

inputs:
  - id: "#compression_level"
    type: ["null", int]
    description: |
      Set compression level, from 0 (uncompressed) to 9 (best)
    inputBinding:
      prefix: "-l"

  - id: "#memory"
    type: ["null", string]
    description: |
      Set maximum memory per thread; suffix K/M/G recognized [768M]
    inputBinding:
      prefix: "-m"

  - id: "#sort_by_name"
    type: ["null", boolean]
    description: "Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates."
    inputBinding:
      prefix: -n

  - id: "#threads"
    type: ["null", int]
    description: "Set number of sorting and compression threads [1]"
    inputBinding:
      prefix: -@

  - id: "#output_name"
    type: string
    description: "Desired output filename."
    inputBinding:
      position: 2

  - id: "#input"
    type: File
    description:
      Input bam file.
    inputBinding:
      position: 1

outputs:
  - id: "#output_file"
    type: File
    outputBinding:
      glob:
        engine: "cwl:JsonPointer"
        script: "job/output_name"

baseCommand: ["samtools", "sort"]

arguments:
  - "-f"
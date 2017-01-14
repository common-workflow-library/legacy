#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

description: "Invoke 'samtools sort' (samtools 1.19)"

inputs:
  compression_level:
    type: int?
    description: |
      Set the desired compression level for the final output file, ranging from
      0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best
      compression but slowest to write), similarly to gzip(1)'s compression
      level setting.

      If -l is not used, the default compression level will apply.
    inputBinding:
      prefix: -l

  memory:
    type: int?
    description: |
      Approximately the maximum required memory per thread, specified  in
      bytes.
    inputBinding:
      prefix: -m

  sort_by_name:
    type: boolean?
    description: "Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates."
    inputBinding:
      prefix: -n

  output_name:
    type: string
    description: "Desired output filename."
    inputBinding:
      position: 2

  input:
    type: File
    description:
      Input bam file.
    inputBinding:
      position: 1

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs['output_name'])

baseCommand: ["samtools", "sort"]

arguments:
  - "-f"


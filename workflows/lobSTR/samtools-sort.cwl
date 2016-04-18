#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"

class: CommandLineTool

description: "Invoke 'samtools sort' (samtools 1.19)"

inputs:
  - id: compression_level
    type: ["null", int]
    description: |
      Set the desired compression level for the final output file, ranging from
      0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best
      compression but slowest to write), similarly to gzip(1)'s compression
      level setting.

      If -l is not used, the default compression level will apply.
    inputBinding:
      prefix: "-l"

  - id: memory
    type: ["null", int]
    description: |
      Approximately the maximum required memory per thread, specified  in
      bytes.
    inputBinding:
      prefix: "-m"

  - id: sort_by_name
    type: ["null", boolean]
    description: "Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates."
    inputBinding:
      prefix: -n

  - id: output_name
    type: string
    description: "Desired output filename."
    inputBinding:
      position: 2

  - id: input
    type: File
    description:
      Input bam file.
    inputBinding:
      position: 1

outputs:
  - id: output_file
    type: File
    outputBinding:
      glob: $(inputs['output_name'])

baseCommand: ["samtools", "sort"]

arguments:
  - "-f"

#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool


inputs:
  output_name:
    type: string
    inputBinding:
      position: 2

    doc: Desired output filename.
  compression_level:
    type: int?
    inputBinding:
      prefix: -l
    doc: |
      Set the desired compression level for the final output file, ranging from
      0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best
      compression but slowest to write), similarly to gzip(1)'s compression
      level setting.

      If -l is not used, the default compression level will apply.
  input:
    type: File
    inputBinding:
      position: 1

    doc: Input bam file.
  sort_by_name:

    type: boolean?
    inputBinding:
      prefix: -n

    doc: Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
  memory:

    type: int?
    inputBinding:
      prefix: -m
    doc: |
      Approximately the maximum required memory per thread, specified  in
      bytes.
outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs['output_name'])

baseCommand: [samtools, sort]

arguments:
- -f
doc: Invoke 'samtools sort' (samtools 1.19)


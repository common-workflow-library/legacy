#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool


requirements:
- class: InitialWorkDirRequirement
  listing:
  - entry: $(inputs['input'])
    entryname: indexed.bam
inputs:
  input:
    type: File
    doc: Input bam file.
outputs:
  bam_with_bai:
    type: File
    outputBinding:
      glob: indexed.bam
    secondaryFiles:
    - .bai

baseCommand: [samtools, index]

arguments:
- indexed.bam
doc: Invoke 'samtools index' to create a 'BAI' index (samtools 1.19)


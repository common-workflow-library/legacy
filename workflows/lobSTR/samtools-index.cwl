#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

description: "Invoke 'samtools index' to create a 'BAI' index (samtools 1.19)"

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

inputs:
  input:
    type: File
    description:
      Input bam file.
    inputBinding:
      valueFrom: $(self.basename)

outputs:
  bam_with_bai:
    type: File
    outputBinding:
      glob: $(inputs.input.basename)
    secondaryFiles:
      - .bai

baseCommand: [samtools, index]


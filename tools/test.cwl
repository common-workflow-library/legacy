#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: ExpressionTool

requirements:
 - class: InlineJavascriptRequirement

inputs:
  bam: File
  bam-index: File

outputs:
  bam-with-index: File

expression: >
  ${
  inputs.bam.secondaryFiles = [ inputs["bam-index"] ];
  return { "bam-with-index": inputs.bam };
  }

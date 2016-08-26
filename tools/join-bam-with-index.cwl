#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: ExpressionTool

doc: Produces a BAM with with its index as a secondaryFile from two distinct
 inputs.

requirements:
 - class: InlineJavascriptRequirement

inputs:
  bam:
   type: File
   format: http://edamontology.org/format_2572
  bam-index: File

outputs:
  bam-with-index: File

expression: >
  ${
  inputs.bam.secondaryFiles = [ inputs["bam-index"] ];
  return { "bam-with-index": inputs.bam }; 
  }

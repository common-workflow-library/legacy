#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/snpeff:4.3

requirements:
  - class: InlineJavascriptRequirement

inputs:

  data_dir:
    type: Directory
    inputBinding:
      prefix: "-dataDir"
      position: 1

  no_stats:
    type: boolean?
    inputBinding:
      prefix: "-noStats"
      position: 2

  csvStats:
    type: boolean?
    inputBinding:
      prefix: "-csvStats"
      position: 3

  output_format:
    type:
      type: enum
      symbols: [ vcf, gatk, bed, bedAnn ]
    default: vcf
    inputBinding:
      prefix: -o
      position: 4

  nodownload:
    type: boolean?
    inputBinding:
      prefix: -nodownload
      position: 5

  verbose:
    type: boolean?
    inputBinding:
      prefix: -v
      position: 6

  genome:
    type: string
    inputBinding:
      position: 7

  variant_calling_file:
    type: File
    inputBinding:
      position: 8

stdout: $(inputs.variant_calling_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.ann.vcf')

outputs:
  annotated_vcf:
    type: stdout

  summary_html:
    type: File?
    outputBinding:
      glob: "snpEff_summary.html"

  summary_txt:
    type: File?
    outputBinding:
      glob: "snpEff_genes.txt"

baseCommand: [ snpEff ]
arguments: [ "ann", "-stats", "snpEff_summary.html" ]

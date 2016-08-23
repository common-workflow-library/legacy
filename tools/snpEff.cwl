#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/snpeff:4.3

requirements:
  - class: InlineJavascriptRequirement

#stdout: $(inputs.inputfile.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.ann.vcf')

inputs:

  genome:
    type: string?
    #default: "GRCh37.75"
    inputBinding:
      position: 1

  variant_calling_file:
    type: File
    format: "http://edamontology.org/format_3016"
    inputBinding:
      position: 2

  genome_dir:
    type: Directory
    inputBinding:
      position: 2

  no_stats:
    type: boolean?
    inputBinding:
      prefix: "-noStats"

  csvStats:
    type: boolean?
    inputBinding:
      prefix: "-csvStats"

  output_format:
    type:
      type: enum
      symbols: [ vcf, gatk, bed, bedAnn ]
    inputBinding:
      prefix: -o

  nodownload:
    type: boolean?
    inputBinding:
      prefix: -nodownload

  verbose:
    type: boolean?
    inputBinding:
      prefix: -v

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
arguments: [ "-stats", "snpEff_summary.html" ]

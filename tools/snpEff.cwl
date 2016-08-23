#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

stdout: $(inputs.inputfile.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.ann.vcf')

inputs:

  genome:
    type: string
    default: "GRCh37.75"
    inputBinding:
      position: 1

  variant_calling_file:
    type: File
    format: "http://edamontology.org/format_3016"
    inputBinding:
      position: 2

  stats_filename:
    type: string?
    default: "snpEff_summary.html"
    inputBinding:
      prefix: "-stats"

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
      glob: $(input.stats_filename)

  summary_txt:
    type: File?
    outputBinding:
      glob: "snpEff_genes.txt"

baseCommand: [ "java", "-Xmx4g" ]

arguments: [ "-jar", "snpEff.jar" ]


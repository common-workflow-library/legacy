#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool


doc: |
      GATK-VariantsFiltration.cwl is developed for CWL consortium

      Usage:
      ```
      java -Xmx8G \
            -jar gatk.jar
            -T VariantFiltration \
            -R reference.fa \
            -V raw_indels.vcf \
            --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
            --filterName "my_indel_filter" \
            -o filtered_indels.vcf
      ```

hints:
- $import: envvar-global.yml
- $import: GATK-docker.yml

inputs:
  reference:
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    inputBinding:
      position: 5
      prefix: -R
    doc: reference genome

  indels_vcf:
    type: File
    inputBinding:
      position: 6
      prefix: -V
    doc: input vcf File raw indel variant calls from haplotype caller

  filter_expression:
    type: string
    default: "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
    inputBinding:
      position: 7
      prefix: --filterExpression

  filter_name:
    type: string
    default: "indel_filter"
    inputBinding:
      position: 8
      prefix: --filterName

  output_filename:
    type: string
    default: "filtered.indels.vcf"
    inputBinding:
      position: 9
      prefix: -o

  java_arg:
    type: string
    default: -Xmx8g
    inputBinding:
      position: 1

  java_tmp:
    type: string
    default: -Djava.io.tmpdir=/tmp
    inputBinding:
      position: 2

outputs:
    output_File:
      type: File
      outputBinding:
        glob: $(inputs.output_filename)
      doc: the filtered vcf file

arguments:
- valueFrom: /usr/local/bin/GenomeAnalysisTK.jar
  position: 3
  prefix: -jar

- valueFrom: VariantFiltration
  position: 4
  prefix: -T

baseCommand: [java]

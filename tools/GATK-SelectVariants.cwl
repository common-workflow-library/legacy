#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool


doc: |
      GATK-SelectVariants.cwl is developed for CWL consortium

      Usage:
      ```
      java -Xmx8G \
            -jar gatk.jar
            -T SelectVariants \
            -R reference.fa \
            -V raw_HC_variants.vcf \
            -selectType INDEL \
            -o raw_indels.vcf
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

  raw_vcf:
    type: File
    inputBinding:
      position: 6
      prefix: -V
    doc: input vcf File raw variant calls from haplotype caller

  select_type:
    type: string
    default: INDEL
    inputBinding:
      position: 7
      prefix: -selectType

  output_filename:
    type: string
    default: "selected.indels.vcf"
    inputBinding:
      position: 8
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
      doc: the selected variants as a vcf file

arguments:
- valueFrom: /usr/local/bin/GenomeAnalysisTK.jar
  position: 3
  prefix: -jar

- valueFrom: SelectVariants
  position: 4
  prefix: -T

baseCommand: [java]

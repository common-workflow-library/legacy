#!/usr/bin/env cwl-runner

# TODO: will need to implement some reasonable (and adjustable) filters

class: Workflow
cwlVersion: v1.0

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  reference:
    type: File
    doc: reference human genome file

  indel_mode:
    type: string
    default: 'INDEL'

  snp_mode:
    type: string
    default: 'SNP'

  snpf_genome:
    type: string

  snpf_nodownload:
    type: boolean

  snpf_data_dir:
    type: Directory

  resource_mills:
    type: File

  haplotest_vcf:
    type: File

  resource_hapmap:
    type: File

  resource_omni:
    type: File

  resource_dbsnp:
    type: File

outputs:

  annotated_indels:
    type: File
    outputSource: snpeff_indels/annotated_vcf

steps:

  snpeff_indels:
    run: ../../tools/snpEff.cwl
    in:
      genome: snpf_genome
      variant_calling_file: haplotest_vcf
      nodownload: snpf_nodownload
      data_dir: snpf_data_dir
    out: [ annotated_vcf ]

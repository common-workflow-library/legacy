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

  select_output_filename:
    type: string
    default: 'selected.indels.vcf'

  filter_expression:
    type: string
    default: "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

  filter_output_filename:
    type: string
    default: 'filtered.indels.vcf'

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

  select_indels:
    run: ../../tools/GATK-SelectVariants.cwl
    in:
      raw_vcf: haplotest_vcf
      reference: reference
      select_type: indel_mode
      output_filename: select_output_filename
    out: [output_File]

  filter_indels:
    run: ../../tools/GATK-VariantFiltration.cwl
    in:
      indels_vcf: select_indels/output_File
      reference: reference
      filter_expression: filter_expression
      output_filename: filter_output_filename
    out: [ output_File ]

  snpeff_indels:
    run: ../../tools/snpEff.cwl
    in:
      genome: snpf_genome
      variant_calling_file: filter_indels/output_File
      nodownload: snpf_nodownload
      data_dir: snpf_data_dir
    out: [ annotated_vcf ]

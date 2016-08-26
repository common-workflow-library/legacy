#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  reference:
    type: File
    doc: reference human genome file
 
  known_ref_db:
    type: File[]?
    doc: array of known variant files for realign target creator

  w_inputBam_HaplotypeCaller:
    type: File
    doc: bam file produced after printReads

  gatk_threads:
    type: int
    doc: number of threads

  outputFileName_HaplotypeCaller:
    type: string
    doc: name of Haplotype caller command output file

  dbsnp:
    type: File
    doc: vcf file containing SNP variations used for Haplotype caller
  
  snpf_genome:
    type: string
  
  snpf_nodownload:
    type: boolean

  snpf_data_dir:
    type: Directory

  indels_resource_mills:
    type: File

  indels_resource_dbsnp:
    type: File

outputs:
  recal_File:
    type: File
    outputSource: vqsr_snps/recal_File
 
steps:

  HaplotypeCaller:
    run: ../../tools/GATK-HaplotypeCaller.cwl
    in:
      outputfile_HaplotypeCaller: outputFileName_HaplotypeCaller
      inputBam_HaplotypeCaller: w_inputBam_HaplotypeCaller
      reference: reference
      dbsnp: dbsnp
      threads: gatk_threads
    out: [ output_HaplotypeCaller ]

  vqsr_snps:
    run: ../../tools/GATK-VariantRecalibrator-SNPs.cwl
    in:
      haplotypecaller_snps_vcf: HaplotypeCaller/output_HaplotypeCaller
      reference: reference
      resource_db: known_ref_db
    out: [tranches_File, recal_File, vqsr_rscript]
  
  vqsr_indels:
    run: ../../tools/GATK-VariantRecalibrator-Indels.cwl
    in:
      haplotypecaller_snps_vcf: HaplotypeCaller/output_HaplotypeCaller
      reference: reference
      resource_mills: indels_resource_mills
      resource_dbsnp: indels_resource_dbsnp
    out: [tranches_File, recal_File, vqsr_rscript]

  apply_recalibration_snps:
    run: ../../tools/GATK-ApplyRecalibration.cwl
    in:
      raw_vcf: HaplotypeCaller/output_HaplotypeCaller
      reference: reference
      recal_file: vqsr_snps/recal_File
      tranches_file: vqsr_snps/tranches_File
    out: [ vqsr_vcf ]

  apply_recalibration_indels:
    run: ../../tools/GATK-ApplyRecalibration.cwl
    in:
      raw_vcf: HaplotypeCaller/output_HaplotypeCaller
      reference: reference
      recal_file: vqsr_indels/recal_File
      tranches_file: vqsr_indels/tranches_File
    out: [ vqsr_vcf ]

  snpeff_snps:
    run: ../../tools/snpEff.cwl
    in:
      genome: snpf_genome
      variant_calling_file: apply_recalibration_snps/vqsr_vcf
      nodownload: snpf_nodownload
      data_dir: snpf_data_dir
    out: [ annotated_vcf ]

  snpeff_indels:
    run: ../../tools/snpEff.cwl
    in:
      genome: snpf_genome
      variant_calling_file: apply_recalibration_indels/vqsr_vcf
      nodownload: snpf_nodownload
      data_dir: snpf_data_dir
    out: [ annotated_vcf ]

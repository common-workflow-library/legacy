#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: Workflow

inputs:
  wp1:
    type: File[]?
    description: list of files containing the first end of paired end reads in fasta or fastq format

  wp2:
    type: File[]?
    description: list of files containing the second end of paired end reads in fasta or fastq format

  w_output_prefix:
    type: string
    description: prefix for output files. will output prefix.aligned.bam and prefix.aligned.stats

  w_reference:
    type: File
    description: "lobSTR's bwa reference files"

  w-rg-sample:
    type: string
    description: Use this in the read group SM tag

  w-rg-lib:
    type: string
    description: Use this in the read group LB tag

  w-strinfo:
    type: File
    description: File containing statistics for each STR.

  w_noise_model:
    type: File
    description: File to read noise model parameters from (.stepmodel)
    secondaryFiles:
      - "^.stuttermodel"

outputs:
  bam:
    type: File
    source: samindex/bam_with_bai

  bam_stats:
    type: File
    source: lobSTR/bam_stats

  vcf:
    type: File
    source: allelotype/vcf

  vcf_stats:
    type: File
    source: allelotype/vcf_stats

hints:
  DockerRequirement:
    dockerLoad: https://workbench.qr1hi.arvadosapi.com/collections/download/qr1hi-4zz18-x2ae13tsx5jqg8d/1nduktd8dpvhdpgsva82lje0i710kgzb6rttks5jldx7s2y7k9/7e0c0ae3bf4e70442f9b8eee816ec23426d9e1169a2925316e5c932745e21613.tar
    dockerImageId: 7e0c0ae3bf4e70442f9b8eee816ec23426d9e1169a2925316e5c932745e21613
    dockerPull: rabix/lobstr

steps:
  lobSTR:
    run: lobSTR-tool.cwl
    in:
      p1: wp1
      p2: wp2
      output_prefix: w_output_prefix
      reference: w_reference
      rg-sample: w-rg-sample
      rg-lib: w-rg-lib
    out: [bam, bam_stats]

  samsort:
    run: samtools-sort.cwl
    in:
      input: lobSTR/bam
      output_name: {default: "aligned.sorted.bam" }
    out: [output_file]

  samindex:
    run: samtools-index.cwl
    in:
      input: samsort/output_file
    out: [bam_with_bai]

  allelotype:
    run: allelotype.cwl
    in:
      bam: samindex/bam_with_bai
      reference: w_reference
      output_prefix: w_output_prefix
      noise_model: w_noise_model
      strinfo: w-strinfo
    out: [vcf, vcf_stats]


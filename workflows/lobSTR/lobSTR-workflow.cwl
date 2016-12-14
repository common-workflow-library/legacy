#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: Workflow

inputs:
  p1:
    type: File[]?
    description: list of files containing the first end of paired end reads in fasta or fastq format

  p2:
    type: File[]?
    description: list of files containing the second end of paired end reads in fasta or fastq format

  output_prefix:
    type: string
    description: prefix for output files. will output prefix.aligned.bam and prefix.aligned.stats

  reference:
    type: File
    description: "lobSTR's bwa reference files"

  rg-sample:
    type: string
    description: Use this in the read group SM tag

  rg-lib:
    type: string
    description: Use this in the read group LB tag

  strinfo:
    type: File
    description: File containing statistics for each STR.

  noise_model:
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
    # dockerLoad: https://workbench.qr1hi.arvadosapi.com/collections/download/qr1hi-4zz18-x2ae13tsx5jqg8d/1nduktd8dpvhdpgsva82lje0i710kgzb6rttks5jldx7s2y7k9/7e0c0ae3bf4e70442f9b8eee816ec23426d9e1169a2925316e5c932745e21613.tar
    # dockerImageId: 7e0c0ae3bf4e70442f9b8eee816ec23426d9e1169a2925316e5c932745e21613
    dockerPull: rabix/lobstr

steps:
  lobSTR:
    run: lobSTR-tool.cwl
    in:
      p1: p1
      p2: p2
      output_prefix: output_prefix
      reference: reference
      rg-sample: rg-sample
      rg-lib: rg-lib
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
      reference: reference
      output_prefix: output_prefix
      noise_model: noise_model
      strinfo: strinfo
    out: [vcf, vcf_stats]


#!/usr/bin/env cwl-runner

class: Workflow

inputs:
  - id: "#p1"
    type:
      - "null"
      - type: array
        items: File
    description: list of files containing the first end of paired end reads in fasta or fastq format

  - id: "#p2"
    type:
      - "null"
      - type: array
        items: File
    description: list of files containing the second end of paired end reads in fasta or fastq format

  - id: "#output_prefix"
    type: string
    description: prefix for output files. will output prefix.aligned.bam and prefix.aligned.stats

  - id: "#reference"
    type: File
    description: "lobSTR's bwa reference files"

  - id: "#rg-sample"
    type: string
    description: Use this in the read group SM tag

  - id: "#rg-lib"
    type: string
    description: Use this in the read group LB tag

  - id: "#strinfo"
    type: File
    description: |
      File containing statistics for each STR.

  - id: "#noise_model"
    type: File
    description: |
      File to read noise model parameters from (.stepmodel)
    inputBinding:
      secondaryFiles:
        - "^.stuttermodel"

outputs:
  - id: "#bam"
    type: File
    source: "#samindex.bam_with_bai"

  - id: "#bam_stats"
    type: File
    source: "#lobSTR.bam_stats"

  - id: "#vcf"
    type: File
    source: "#allelotype.vcf"

  - id: "#vcf_stats"
    type: File
    source: "#allelotype.vcf_stats"

hints:
  - class: DockerRequirement
    dockerLoad: https://workbench.qr1hi.arvadosapi.com/collections/download/qr1hi-4zz18-zjbmbgf5felkl4q/8g3zdu8nwv9emx21jlw1syovy4tfuzwl6ym26j51307rw8vkv/b2d7226fe5130f01cc3120bfb69a6f9c5afbfec1f175ddb24e8200abf4c7e8ed.tar
    dockerImageId: b2d7226fe5130f01cc3120bfb69a6f9c5afbfec1f175ddb24e8200abf4c7e8ed

steps:
  - id: "#lobSTR"
    run: { import: lobSTR-tool.cwl }
    inputs:
      - { id: "#lobSTR.p1", source: "#p1" }
      - { id: "#lobSTR.p2", source: "#p2" }
      - { id: "#lobSTR.output_prefix", source: "#output_prefix" }
      - { id: "#lobSTR.reference", source: "#reference" }
      - { id: "#lobSTR.rg-sample", source: "#rg-sample" }
      - { id: "#lobSTR.rg-lib", source: "#rg-lib" }
    outputs:
      - { id: "#lobSTR.bam" }
      - { id: "#lobSTR.bam_stats" }

  - id: "#samsort"
    run: { import: samtools-sort.cwl }
    inputs:
      - { id: "#samsort.input", source: "#lobSTR.bam" }
      - { id: "#samsort.output_name", default: "aligned.sorted.bam" }
    outputs:
      - { id: "#samsort.output_file" }

  - id: "#samindex"
    run: { import: samtools-index.cwl }
    inputs:
      - { id: "#samindex.input", source: "#samsort.output_file" }
    outputs:
      - { id: "#samindex.bam_with_bai" }

  - id: "#allelotype"
    run: { import: allelotype.cwl }
    inputs:
      - { id: "#allelotype.bam", source: "#samindex.bam_with_bai" }
      - { id: "#allelotype.reference", source: "#reference" }
      - { id: "#allelotype.output_prefix", source: "#output_prefix" }
      - { id: "#allelotype.noise_model", source: "#noise_model" }
      - { id: "#allelotype.strinfo", source: "#strinfo" }
    outputs:
      - { id: "#allelotype.vcf" }
      - { id: "#allelotype.vcf_stats" }

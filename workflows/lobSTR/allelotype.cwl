#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

description: Run lobSTR allelotype classifier.

requirements:
 InlineJavascriptRequirement: {}

inputs:
  bam:
    type: File
    description: |
      BAM file to analyze. Should have a unique read group and be sorted and indexed.
    inputBinding:
      prefix: "--bam"
    secondaryFiles:
      - ".bai"

  output_prefix:
    type: string
    description: "Prefix for output files. will output prefix.vcf and prefix.genotypes.tab"
    inputBinding:
      prefix: "--out"

  noise_model:
    type: File
    description: |
      File to read noise model parameters from (.stepmodel)
    inputBinding:
      prefix: "--noise_model"
      valueFrom: |
          $( {"path": self.dirname + "/" + self.nameroot, "class": "File"})
    secondaryFiles:
      - "^.stuttermodel"

  strinfo:
    type: File
    description: |
      File containing statistics for each STR.
    inputBinding:
      prefix: "--strinfo"

  reference:
    type: File
    description: "lobSTR's bwa reference files"
    inputBinding:
      prefix: "--index-prefix"
      valueFrom: |
          $({"path": self.path.replace(/ref\.fasta$/, ""), "class": "File"})

    secondaryFiles:
      - ".amb"
      - ".ann"
      - ".bwt"
      - ".pac"
      - ".rbwt"
      - ".rpac"
      - ".rsa"
      - $(self.path.replace(/ref$/, "chromsizes.tab"))
      - $(self.path.replace(/ref$/, "mergedref.bed"))
      - $(self.path)_map.tab

outputs:
  vcf:
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.vcf')
  vcf_stats:
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.allelotype.stats')

baseCommand: ["allelotype", "--command", "classify"]

arguments:
  - "--noweb"

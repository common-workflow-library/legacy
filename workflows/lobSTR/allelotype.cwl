#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"

class: CommandLineTool

description: Run lobSTR allelotype classifier.

requirements:
 - class: InlineJavascriptRequirement

inputs:
  - id: bam
    type: File
    description: |
      BAM file to analyze. Should have a unique read group and be sorted and indexed.
    inputBinding:
      prefix: "--bam"
    secondaryFiles:
      - ".bai"

  - id: output_prefix
    type: string
    description: "Prefix for output files. will output prefix.vcf and prefix.genotypes.tab"
    inputBinding:
      prefix: "--out"

  - id: noise_model
    type: File
    description: |
      File to read noise model parameters from (.stepmodel)
    inputBinding:
      prefix: "--noise_model"
      valueFrom: |
          ${ return {"path": self.path.match(/(.*)\.stepmodel/)[1], "class": "File"}; }
    secondaryFiles:
      - "^.stuttermodel"

  - id: strinfo
    type: File
    description: |
      File containing statistics for each STR.
    inputBinding:
      prefix: "--strinfo"

  - id: reference
    type: File
    description: "lobSTR's bwa reference files"
    inputBinding:
      prefix: "--index-prefix"
      valueFrom: |
          ${ return {"path": self.path.match(/(.*)ref\.fasta/)[1], "class": "File"}; }

    secondaryFiles:
      - ".amb"
      - ".ann"
      - ".bwt"
      - ".pac"
      - ".rbwt"
      - ".rpac"
      - ".rsa"
      - ${return self.path.replace(/(.*)ref\.fasta/, "$1chromsizes.tab")}
      - ${return self.path.replace(/(.*)ref\.fasta/, "$1mergedref.bed")}
      - ${return self.path.replace(/(.*)ref\.fasta/, "$1ref_map.tab")}

outputs:
  - id: vcf
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.vcf')
  - id: "#vcf_stats"
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.allelotype.stats')

baseCommand: ["allelotype", "--command", "classify"]

arguments:
  - "--noweb"

#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3.dev2"

class: CommandLineTool

description: |
  Usage: vcftobedbpe -i <in.vcf> -o [out.bedpe]

requirements:
  - "@import": envvar-global.cwl

inputs:
  - id: "#input"
    type: File
    description: |
      "Input vcf file."
    inputBinding:
      prefix: "-i"


outputs:
  - id: "#output"
    type: File
    description: "The bedpe file"
    outputBinding:
      prefix: "-o"

baseCommand: ["samtools"]

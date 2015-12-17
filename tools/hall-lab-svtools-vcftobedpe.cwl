#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3.dev2"

class: CommandLineTool

description: |
  Usage: vcftobedbpe -i <in.vcf> -o [out.bedpe]

requirements:
  - "@import": envvar-global.cwl

hints:
  - class: DockerRequirement
    dockerPull: biocrusoe/hall-lab-svtools

inputs:
  - id: "#input"
    type: File
    description: |
      "Input vcf file."
    streamable: true
    inputBinding:
      prefix: "-i"

stdout:
  "output.bedpe"

outputs:
  - id: "#bedpe"
    type: File
    description: "The bedpe file"
    streamable: true
    outputBinding:
      glob: "output.bedpe"

baseCommand: ["vcftobedpe"]

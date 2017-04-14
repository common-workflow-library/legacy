#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
 DockerRequirement:
  dockerPull: insilicodb/kallisto

inputs:
 fastqs:
   type: File[]
   format: http://edamontology.org/format_1930 # FASTA
   inputBinding: {}

 index:
   type: File
   inputBinding:
     prefix: "--index"

baseCommand: [ kallisto, quant ]

arguments: [ "--output-dir", out ]

outputs:
 quantification:
  type: File
  outputBinding:
   glob: out/abundance.h5
   

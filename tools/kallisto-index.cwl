#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
 DockerRequirement:
  dockerPull: insilicodb/kallisto

inputs:
 fasta-files:
   type: File[]
   format: http://edamontology.org/format_1929 # FASTA
   inputBinding: {}

baseCommand: [ kallisto, index ]

arguments: [ --index, index.idx ]

outputs:
 index:
  type: File
  outputBinding:
   glob: $(inputs.index_name)
   

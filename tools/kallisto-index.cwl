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
   doc: 
   inputBinding:
    position: 100

 index_name:
   type: string?
   default: index.idx
   inputBinding:
     prefix: "--index"

baseCommand: [ kallisto, index ]

outputs:
 index:
  type: File
  outputBinding:
   glob: $(inputs.index_name)
   

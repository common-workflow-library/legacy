#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

requirements:
  - $import: envvar-global.yml                                                      
  - $import: bwa-docker.yml

  - class: InlineJavascriptRequirement

inputs:
  prefix:
  type: File
  inputBinding:
    position: 4

  input:
    type: File
    inputBinding: 5

  output_filename:
    type: string
    inputBinding:
      position: 1
      prefix: "-f"

  threads:
    type: int?
    inputBinding:
      postion: 1
      prefix: "-t"

  n:
    type: [ "null",int,float ]
    inputBinding:
      position: 1
      prefix: "-n"

  o:
    type: int?
    inputBinding:
      position: 1
      prefix: "-o"

  e: 
    type: int?
    inputBinding:
      prefix: "-e"
  
  i:
    type: int?
    inputBinding:
      prefix: "-i"

  d:
    type: int?
    inputBinding:
      prefix: "-d"

  l:
    type: int?
    inputBinding:
      prefix: "-l"

  k:
    type: int?
    inputBinding:
      prefix: "-k"

  m:
    type: int?
    inputBinding:
      prefix: "-m"

  I: 
    type: boolean?
    inputBinding:
      prefix: "-I"

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

baseCommand: [bwa, aln]


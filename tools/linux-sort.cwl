#!/usr/bin/env cwl-runner

class: CommandLineTool

description: |
  Usage: sort [OPTION]... [FILE]...
    or:  sort [OPTION]... --files0-from=F
  Write sorted concatenation of all FILE(s) to standard output.

requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: linux-sort-docker.cwl

inputs:
  - id: "#input"
    type: File
    inputBinding:
      position: 4

  - id: "#key"
    type: 
      type: array
      items: string 
    description: |
      -k, --key=POS1[,POS2]
      start a key at POS1, end it at POS2 (origin 1)

outputs:
  - id: "#sorted"
    type: File
    description: "The sorted file"
    outputBinding:
      glob:
        engine: node-engine.cwl
        script: $job['input'].path.split('/').slice(-1)[0] + '.sorted'

stdout: 
  engine: node-engine.cwl
  script: $job['input'].path.split('/').slice(-1)[0] + '.sorted'

baseCommand: ["sort"]

arguments:
  - valueFrom:
      engine: node-engine.cwl
      script: $job['key'].map(function(i) {return "-k"+i;})


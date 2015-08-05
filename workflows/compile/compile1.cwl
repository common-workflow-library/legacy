#!/usr/bin/env cwl-runner


- id: "#compile"
  class: CommandLineTool
  inputs:
    - id: "#src"
      type: File
      inputBinding: {}
    - id: "#object"
      type: string
      inputBinding:
          prefix: "-o"
  outputs:
    - id: "#compiled"
      type: File
      outputBinding:
          glob:
             engine: cwl:JsonPointer
             script: /job/object
  baseCommand: gcc
  arguments:
     - "-c"
     - "-Wall"

- id: "#link"
  class: CommandLineTool
  inputs:
    - id: "#objects"
      type:  { type: array, items: File }
      inputBinding: {}
    - id: "#output"
      type: string
      inputBinding:
          prefix: "-o"
  outputs:
    - id: "#executable"
      type: File
      outputBinding:
          glob:
             engine: cwl:JsonPointer
             script: /job/output
  baseCommand: gcc


- id: "#main"
  class: Workflow
  inputs: []
  outputs:
    - id: "#main.output"
      type: File
      source: "#compile.compiled"
  steps :
    - id: "#compilesources-src1"
      run: {import: "#compile"}
      inputs:
         - { id: "#compile.src" , default: "source1.c" }
         - { id: "#compile.object" , default: "source1.o" }
      outputs:
        - { id: "#compile.compiled" }
    - id: "#compilesources-src2"
      run: {import: "#compile"}
      inputs:
         - { id: "#compile.src" , default: "source2.c" }
         - { id: "#compile.object" , default: "source2.o" }
      outputs:
         - { id: "#compile.compiled" }


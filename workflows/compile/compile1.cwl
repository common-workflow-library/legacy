#!/usr/bin/env cwl-runner

- id: "#compile"
  cwlVersion: "cwl:draft-3"
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
  cwlVersion: "cwl:draft-3"
  inputs:
    - id: "#objects"
      type:  { type: array, items: File }
      inputBinding:
        position: 2
    - id: "#output"
      type: string
      inputBinding:
          position: 1
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
  cwlVersion: "cwl:draft-3"
  requirements:
    - class: EnvVarRequirement
      envDef:
        - envName: PATH
          envValue: /usr/bin
  inputs: []
  outputs:
    - id: "#main.output"
      type: File
      source: "#linkobj.executable"
  steps :
    - id: "#compilesources-src1"
      run: {import: "#compile"}
      inputs:
         - { id: "#compilesources-src1.src" , default: {class: File, path: "source1.c" } }
         - { id: "#compilesources-src1.object" , default: "source1.o" }
      outputs:
        - { id: "#compilesources-src1.compiled" }

    - id: "#compilesources-src2"
      run: {import: "#compile"}
      inputs:
         - { id: "#compilesources-src2.src" , default: {class: File, path: "source2.c" } }
         - { id: "#compilesources-src2.object" , default: "source2.o" }
      outputs:
         - { id: "#compilesources-src2.compiled" }

    - id: "#linkobj"
      run: {import: "#link"}
      inputs:
         - { id: "#linkobj.objects" , source: ["#compilesources-src1.compiled", "#compilesources-src2.compiled"] }
         - { id: "#linkobj.output" , default: "a.out" }
      outputs:
         - { id: "#linkobj.executable" }


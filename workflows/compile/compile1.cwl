#!/usr/bin/env cwl-runner

cwlVersion: v1.0
$graph:
- id: compile
  class: CommandLineTool
  inputs:
    src:
      type: File
      inputBinding: {}
    object:
      type: string
      inputBinding:
          prefix: "-o"
  outputs:
    compiled:
      type: File
      outputBinding:
        glob: $(inputs.object)
  baseCommand: gcc
  arguments:
    - "-c"
    - "-Wall"

- id: link
  class: CommandLineTool
  inputs:
    objects:
      type:  File[]
      inputBinding:
        position: 2
    output:
      type: string
      inputBinding:
          position: 1
          prefix: "-o"
  outputs:
    executable:
      type: File
      outputBinding:
          glob: $(inputs.output)
  baseCommand: gcc

- id: main
  class: Workflow
  requirements:
    - class: MultipleInputFeatureRequirement
  inputs: []
  outputs:
    - id: output
      type: File
      outputSource: linkobj/executable
  steps:
    compilesources-src1:
      run: "#compile"
      in:
         src:
           default:
             class: File
             location: source1.c
             secondaryFiles:
               - class: File
                 location: source1.h
         object: { default: "source1.o" }
      out: [compiled]

    compilesources-src2:
      run: "#compile"
      in:
         src: { default: {class: File, location: "source2.c" } }
         object: { default: "source2.o" }
      out: [compiled]

    linkobj:
      run: "#link"
      in:
         objects: [compilesources-src1/compiled, compilesources-src2/compiled]
         output: { default: "a.out" }
      out: [executable]

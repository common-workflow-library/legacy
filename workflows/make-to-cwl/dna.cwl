#!/usr/bin/env cwl-runner

- id: rnaseq
  cwlVersion: "cwl:draft-3"
  class: CommandLineTool
  inputs:
    - id: sequence
      type: string
      inputBinding: {}
  outputs:
    - id: seqout
      type: File
      outputBinding:
        glob: rnaseq
  baseCommand: echo
  stdout: rnaseq

- id: cat
  cwlVersion: "cwl:draft-3"
  class: CommandLineTool
  inputs:
    - id: sequences
      type:
        type: array
        items: File
      inputBinding: {}
    - id: catfilename
      type: string
  outputs:
    - id: catout
      type: File
      outputBinding:
        glob: $(inputs["catfilename"])
  baseCommand: cat
  stdout: $(inputs["catfilename"])

- id: tr
  cwlVersion: "cwl:draft-3"
  class: CommandLineTool
  inputs:
    - id: trinput
      type: File
    - id: from
      type: string
      inputBinding:
        position: 1
    - id: to
      type: string
      inputBinding:
        position: 1
    - id: filename
      type: string
  outputs:
    - id: trout
      type: File
      outputBinding:
        glob: $(inputs.filename)
  baseCommand: tr
  stdin: $(inputs.trinput.path)
  stdout: $(inputs.filename)

- id: main
  cwlVersion: "cwl:draft-3"
  class: Workflow
  inputs:
    - id: rna
      type:
        type: array
        items: string
  outputs:
    - id: outfile
      type: File
      source: "#main/combine_sequences/catout"

  requirements:
    - class: ScatterFeatureRequirement

  steps:
    - id: get_sequences
      run: "#rnaseq"
      scatter: "#main/get_sequences/sequence"
      inputs:
        - { id: sequence, source: "#main/rna" }
      outputs:
        - { id: seqout }

    - id: translate_sequences
      run: "#tr"
      scatter: "#main/translate_sequences/trinput"
      inputs:
        - { id: trinput, source: "#main/get_sequences/seqout" }
        - { id: from, default: "U" }
        - { id: to, default: "T" }
        - { id: "#combine_sequences.filename", default: "dna" }
      outputs:
        - { id: "#translate_sequences.trout" }

    - id: combine_sequences
      run: "#cat"
      inputs:
        - { id: sequences, source: "#translate_sequences.trout" }
        - { id: catfilename, default: "database.dna" }
      outputs:
        - { id: catout }

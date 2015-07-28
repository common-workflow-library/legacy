#!/usr/bin/env cwl-runner

- id: "#rnaseq"
  class: CommandLineTool
  inputs:
    - id: "#sequence"
      type: string
      inputBinding: {}
  outputs:
    - id: "#seqout"
      type: File
      outputBinding:
        glob: rnaseq
  baseCommand: echo
  stdout: rnaseq

- id: "#cat"
  class: CommandLineTool
  inputs:
    - id: "#sequences"
      type:
        type: array
        items: File
      inputBinding: {}
    - id: "#catfilename"
      type: string
  outputs:
    - id: "#catout"
      type: File
      outputBinding:
        glob:
          engine: cwl:JsonPointer
          script: /job/catfilename
  baseCommand: cat
  stdout:
    engine: cwl:JsonPointer
    script: /job/catfilename


- id: "#tr"
  class: CommandLineTool
  inputs:
    - id: "#trinput"
      type: File
    - id: "#from"
      type: string
      inputBinding:
        position: 1
    - id: "#to"
      type: string
      inputBinding:
        position: 1
    - id: "#filename"
      type: string
  outputs:
    - id: "#trout"
      type: File
      outputBinding:
        glob:
          engine: cwl:JsonPointer
          script: /job/filename
  baseCommand: tr
  stdin:
    engine: cwl:JsonPointer
    script: /job/trinput/path
  stdout:
    engine: cwl:JsonPointer
    script: /job/filename


- id: "#main"
  class: Workflow
  inputs:
    - id: "#rna"
      type:
        type: array
        items: string
  outputs:
    - id: "#outfile"
      type: File
      source: "#combine_sequences.catout"

  requirements:
    - class: ScatterFeatureRequirement

  steps:
    - id: "#get_sequences"
      run: {import: "#rnaseq"}
      scatter: "#get_sequences.sequence"
      inputs:
        - { id: "#get_sequences.sequence", source: "#rna" }
      outputs:
        - { id: "#get_sequences.seqout" }

    - id: "#translate_sequences"
      run: {import: "#tr"}
      scatter: "#translate_sequences.trinput"
      inputs:
        - { id: "#translate_sequences.trinput", source: "#get_sequences.seqout" }
        - { id: "#translate_sequences.from", default: "U" }
        - { id: "#translate_sequences.to", default: "T" }
        - { id: "#combine_sequences.filename", default: "dna" }
      outputs:
        - { id: "#translate_sequences.trout" }

    - id: "#combine_sequences"
      run: {import: "#cat"}
      inputs:
        - { id: "#combine_sequences.sequences", source: "#translate_sequences.trout" }
        - { id: "#combine_sequences.catfilename", default: "database.dna" }
      outputs:
        - { id: "#combine_sequences.catout" }

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
$graph:
- id: rnaseq
  class: CommandLineTool
  inputs:
    sequence: string
  outputs:
    seqout: stdout
  arguments: [echo, $(inputs.sequence)]
  stdout: rnaseq

- id: cat
  class: CommandLineTool
  inputs:
    sequences: File[]
    catfilename: string
  outputs:
    catout: stdout
  arguments: [cat, $(inputs.sequences)]
  stdout: $(inputs.catfilename)

- id: tr
  class: CommandLineTool
  inputs:
    trinput: File
    from: string
    to: string
    filename: string
  outputs:
    trout: stdout
  arguments: [tr, $(inputs.from), $(inputs.to)]
  stdin: $(inputs.trinput.path)
  stdout: $(inputs.filename)

- id: main
  class: Workflow
  inputs:
    rna: string[]
  outputs:
    outfile:
      type: File
      outputSource: translate_sequences/trout

  requirements:
    - class: ScatterFeatureRequirement

  steps:
    get_sequences:
      run: "#rnaseq"
      scatter: sequence
      in:
        sequence: rna
      out: [seqout]

    combine_sequences:
      run: "#cat"
      in:
        sequences: get_sequences/seqout
        catfilename: { default: "database.dna" }
      out: [catout]

    translate_sequences:
      run: "#tr"
      in:
        trinput: combine_sequences/catout
        from: { default: "U" }
        to: { default: "T" }
        filename: { default: "database.dna" }
      out: [trout]

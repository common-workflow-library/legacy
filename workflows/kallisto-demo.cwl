#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

# Inspired by
# https://github.com/InSilicoDB/pipeline-kallisto/blob/develop/main.nf

hints:
 DockerRequirement:
  dockerPull: insilicodb/kallisto

inputs:
  transcripts: File[]
  reads: File[]

steps:
  indexing:
    run: ../tools/kallisto-index.cwl
    in:
      fasta-files: transcripts
    out:
     - index
  quantifying:
    run: ../tools/kallisto-quant.cwl
    in:
      index: indexing/index
      fastqs: reads
    out:
      - quantification

outputs:
  quants:
    type: File
    outputSource: quantifying/quantification



#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement

inputs:
  fastqSeqs: File[]

outputs:
  reports:
    type: Directory[]
    outputSource: runFastqc/report

steps:
  runFastqc:
    run: ../../tools/fastqc.cwl
    in:
      fastqFile: fastqSeqs
    scatter: fastqFile
    out: [ report ]

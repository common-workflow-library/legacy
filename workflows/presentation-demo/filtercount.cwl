#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: "cwl:draft-3"

requirements:
  - class: ScatterFeatureRequirement
  - class: DockerRequirement
    dockerPull: "debian:8"

inputs:
  - id: pattern
    type: string
  - id: infile
    type: {type: array, items: File}

outputs:
  - id: outfile
    type: File
    source: "#wc/outfile"

steps:
  - id: grep
    run: grep.cwl
    scatter: "#grep/infile"
    inputs:
      - id: pattern
        source: "#pattern"
      - id: infile
        source: "#infile"
    outputs:
      - id: outfile

  - id: wc
    run: wc.cwl
    inputs:
      - id: infile
        source: "#grep/outfile"
    outputs:
      - id: outfile

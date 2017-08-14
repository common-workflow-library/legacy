#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
  ScatterFeatureRequirement: {}
  DockerRequirement: 
    dockerPull: debian:8

inputs:
  pattern: string
  texts: File[]

outputs:
  outfile:
    type: File
    outputSource: count_matches/counts

steps:
  search_texts:
    run: grep.cwl
    inputs:
      pattern: pattern
      file_to_search: texts
    scatter: file_to_search
    outputs: [ results ]

  count_matches:
    run: wc.cwl
    inputs:
      files: search_texts/results
    outputs: [ counts ]


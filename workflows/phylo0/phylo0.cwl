#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
 - class: ScatterFeatureRequirement

inputs:
  unaligned: File
  phylo_method: string[]

outputs:
  result: File

steps:
  alignSeqs:
    run: emma.cwl
    in:
      unaligned: unaligned
    out: [ alignedseq ]

  computeDist:
    run: fnadist.cwl
    in:
      aligned: alignSeqs/alignedseq
      method: phylo_method
    scatter: method
    out: [ distance ]

  makePhylo:
    run: fneighbor.cwl
    in:
      dist: computeDist/distance
    out: [ result ]

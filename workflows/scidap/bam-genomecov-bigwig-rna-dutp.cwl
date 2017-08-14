#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: Workflow

requirements:
- class: SubworkflowFeatureRequirement

inputs:
- id: input
  type: File

- id: genomeFile
  type: File

- id: scale
  type: float

- id: bigWigP
  type: string

- id: bigWigR
  type: string

outputs:
- id: outfileP
  type: File
  outputSource: '#genomecovP.outfile'

- id: outfileR
  type: File
  outputSource: '#genomecovR.outfile'

steps:
- id: '#genomecovP'
  run: bam-genomecov-bigwig.cwl
  inputs:
  - {id: input, source: '#input'}
  - {id: genomeFile, source: '#genomeFile'}
  - {id: strand, default: +}
  - {id: scale, source: '#scale'}
  - {id: bigWig, source: '#bigWigP'}
  outputs:
  - {id: outfile"}

- id: '#genomecovR'
  run: bam-genomecov-bigwig.cwl
  inputs:
  - {id: input, source: '#input'}
  - {id: genomeFile, source: '#genomeFile'}
  - {id: strand, default: '-'}
  - {id: scale, source: '#scale'}
  - {id: bigWig, source: '#bigWigR'}
  outputs:
  - {id: outfile}


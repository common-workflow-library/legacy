#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"

class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - $import: ../../tools/envvar-global.yml

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
    source: "#genomecovP.outfile"

  - id: outfileR
    type: File
    source: "#genomecovR.outfile"

steps:
  - id: "#genomecovP"
    run: bam-genomecov-bigwig.cwl
    inputs:
      - {id: input, source: "#input"}
      - {id: genomeFile, source: "#genomeFile"}
      - {id: strand, default: "+" }
      - {id: scale, source: "#scale" }
      - {id: bigWig, source: "#bigWigP" }
    outputs:
      - {id: outfile"}

  - id: "#genomecovR"
    run: bam-genomecov-bigwig.cwl
    inputs:
      - {id: input, source: "#input"}
      - {id: genomeFile, source: "#genomeFile"}
      - {id: strand, default: "-" }
      - {id: scale, source: "#scale" }
      - {id: bigWig, source: "#bigWigR" }
    outputs:
      - {id: outfile}

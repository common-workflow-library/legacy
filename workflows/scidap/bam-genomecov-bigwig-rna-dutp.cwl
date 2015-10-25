#!/usr/bin/env cwl-runner

class: Workflow
requirements:
  - class: SubworkflowFeatureRequirement
  - import: node-engine.cwl
  - import: ../../tools/envvar-global.cwl

inputs:
  - id: "#input"
    type: File

  - id: "#genomeFile"
    type: File

  - id: "#scale"
    type: float

  - id: "#bigWigP"
    type: string

  - id: "#bigWigR"
    type: string

outputs:
  - id: "#outfileP"
    type: File
    source: "#genomecovP.outfile"

  - id: "#outfileR"
    type: File
    source: "#genomecovR.outfile"

steps:
  - id: "#genomecovP"
    run: {import: ./bam-genomecov-bigwig.cwl}
    inputs:
      - {id: "#genomecovP.input", source: "#input"}
      - {id: "#genomecovP.genomeFile", source: "#genomeFile"}
      - {id: "#genomecovP.strand", default: "+" }
      - {id: "#genomecovP.scale", source: "#scale" }
      - {id: "#genomecovP.bigWig", source: "#bigWigP" }
    outputs:
      - {id: "#genomecovP.outfile"}

  - id: "#genomecovR"
    run: {import: ./bam-genomecov-bigwig.cwl}
    inputs:
      - {id: "#genomecovR.input", source: "#input"}
      - {id: "#genomecovR.genomeFile", source: "#genomeFile"}
      - {id: "#genomecovR.strand", default: "-" }
      - {id: "#genomecovR.scale", source: "#scale" }
      - {id: "#genomecovR.bigWig", source: "#bigWigR" }
    outputs:
      - {id: "#genomecovR.outfile"}

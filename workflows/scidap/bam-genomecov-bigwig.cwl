#!/usr/bin/env cwl-runner

class: Workflow
requirements:
  - import: node-engine.cwl
  - import: ../../tools/envvar-global.cwl
#  - import: ../../tools/bedtools-genomecov-types.cwl

inputs:
  - id: "#input"
    type: File

  - id: "#genomeFile"
    type: File

  - id: "#scale"
    type: float

  - id: "#bigWig"
    type: string

outputs:
  - id: "#outfile"
    type: File
    source: "#bigwig.bigWigOut"

steps:
  - id: "#genomecov"
    run: {import: ../../tools/bedtools-genomecov.cwl}
    inputs:
      - {id: "#genomecov.input", source: "#input"}
      - {id: "#genomecov.genomeFile", source: "#genomeFile"}
      - {id: "#genomecov.genomecoverageout", default: "genomecov.bed" }
      - {id: "#genomecov.dept", default: '-bg' }
      - {id: "#genomecov.split", default: true }
      - {id: "#genomecov.scale", source: "#scale" }
    outputs:
      - {id: "#genomecov.genomecoverage"}

  - id: "#sort"
    run: {import: ../../tools/linux-sort.cwl}
    inputs:
      - {id: "#sort.input", source: "#genomecov.genomecoverage" }
      - {id: "#sort.key", default: ["1,1","2,2n"] }
    outputs:
      - {id: "#sort.sorted"}

  - id: "#bigwig"
    run: {import: ../../tools/ucsc-bedGraphToBigWig.cwl}
    inputs:
      - {id: "#bigwig.input", source: "#sort.sorted"}
      - {id: "#bigwig.genomeFile", source: "#genomeFile"}
      - {id: "#bigwig.bigWig", source: "#bigWig"}
    outputs:
      - {id: "#bigwig.bigWigOut"}

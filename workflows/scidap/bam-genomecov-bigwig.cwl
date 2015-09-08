#!/usr/bin/env cwl-runner

class: Workflow
inputs:
  - id: "#type"
    type: string
  - id: "#input"
    type: File
  - id: "#genomeFile"
    type: File

outputs:
  - id: "#outfile"
    type: File
    #source: "#bigwig.fileout"
    source: "#genomecov.genomecoverage"

steps:
  - id: "#genomecov"
    run: {import: ../../tools/bedtools-genomecov.cwl}
    inputs:
      - {id: "#genomecov.input", source: "#input"}
      - {id: "#genomecov.genomeFile", source: "#genomeFile"}
      - {id: "#genomecov.genomecoverageout", default: "genomecov.bed" }
      - {id: "#genomecov.dept", default: "{ '#dept': '-bg' }"  }
    outputs:
      - {id: "#genomecov.genomecoverage"}


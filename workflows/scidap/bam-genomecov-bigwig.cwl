#!/usr/bin/env cwl-runner

class: Workflow
requirements:
  - import: ../../engines/node-engine.cwl
  - import: ../../tools/envvar-global.cwl
  - import: ../../tools/bedtools-genomecov-types.cwl

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
      - {id: "#genomecov.dept", type: '#depts' ,default: {'#dept': '-bg'} }
    outputs:
      - {id: "#genomecov.genomecoverage"}


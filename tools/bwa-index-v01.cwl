#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool


requirements:
- $import: envvar-global.yml
- $import: envvar-global.yml
- class: InlineJavascriptRequirement

inputs:
  input:
    type: File
    inputBinding: {}

  a:
    type: string?
    inputBinding:
      prefix: "-a"

  p:
    type: string?
    inputBinding:
      prefix: "-p"

  six:
    type: boolean?
    inputBinding:
      prefix: "-6"

outputs:
  output:
    type:
      type: array
      items: File
    outputBinding:
      glob:
        - ${
              if (inputs.p) {
                return inputs.p + ".amb"
              } else {
                if (inputs._6 == true) {
                  return inputs.input.path + ".64.amb"
                } else {
                  return inputs.input.path + ".amb"
                }
              }
            }
          - ${
              if (inputs.p) {
                return inputs.p + ".ann"
              } else {
                if (inputs._6 == true) {
                  return inputs.input.path + ".64.ann"
                } else {
                  return inputs.input.path + ".ann"
                }
              }
            }
          - ${
              if (inputs.p) {
                return inputs.p + ".bwt"
              } else {
                if (inputs._6 == true) {
                  return inputs.input.path + ".64.bwt"
                } else {
                  return inputs.input.path + ".bwt"
                }
              }
            }
          - ${
              if (inputs.p) {
                return inputs.p + ".pac"
              } else {
                if (inputs._6 == true) {
                  return inputs.input.path + ".64.pac"
                } else {
                  return inputs.input.path + ".pac"
                }
              }
            }
          - ${
              if (inputs.p) {
                return inputs.p + ".sa"
              } else {
                if (inputs._6 == true) {
                  return inputs.input.path + ".64.sa"
                } else {
                  return inputs.input.path + ".sa"
                }
              }
            }

baseCommand:
  - bwa
  - index

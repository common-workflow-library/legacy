#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- class: InlineJavascriptRequirement

#TODO: Enable after this issue is fixed: https://github.com/common-workflow-language/cwltool/issues/80
#hints:
#  - $import: bwa-docker.yml

inputs:
  a:
    type: string?
    inputBinding:
      position: 2
      prefix: -a
    doc: |
      BWT construction algorithm: bwtsw or is (Default: auto)
  input:
    type: File
    inputBinding:
      position: 3

  b:
    type: int?
    inputBinding:
      position: 2
      prefix: -b
    doc: |
      Block size for the bwtsw algorithm (effective with -a bwtsw) (Default: 10000000)
  _6:
    type: boolean?
    inputBinding:
      position: 2
      prefix: '-6'
    doc: |
      Index files named as <in.fasta>.64.* instead of <in.fasta>.*
  p:
    type: string?
    inputBinding:
      position: 2
      prefix: -p
    doc: |
      Prefix of the index (Default: same as fasta name)
outputs:
  output:
    type: {type: array, items: File}
    outputBinding:
      glob:
      - |
        ${
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
      - |
        ${
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
      - |
        ${
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
      - |
        ${
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
      - |
        ${
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

doc: |
  Usage:   bwa index [options] <in.fasta>

  Options: -a STR    BWT construction algorithm: bwtsw or is [auto]
           -p STR    prefix of the index [same as fasta name]
           -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
           -6        index files named as <in.fasta>.64.* instead of <in.fasta>.*

  Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
           `-a div' do not work not for long genomes.


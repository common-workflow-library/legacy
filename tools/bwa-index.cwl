#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

requirements:
  - $import: envvar-global.yml
  - class: InlineJavascriptRequirement

#TODO: Enable after this issue is fixed: https://github.com/common-workflow-language/cwltool/issues/80
#hints:
#  - $import: bwa-docker.yml

inputs:
  - id: "input"
    type: File
    inputBinding:
      position: 3

  - id: "a"
    type: ["null", string]
    description: |
      BWT construction algorithm: bwtsw or is (Default: auto)
    inputBinding:
      position: 2
      prefix: "-a"

  - id: "p"
    type: ["null", string]
    description: |
      Prefix of the index (Default: same as fasta name)
    inputBinding:
      position: 2
      prefix: "-p"

  - id: "b"
    type: ["null", int]
    description: |
      Block size for the bwtsw algorithm (effective with -a bwtsw) (Default: 10000000)
    inputBinding:
      position: 2
      prefix: "-b"

  - id: "_6"
    type: ["null", boolean]
    description: |
      Index files named as <in.fasta>.64.* instead of <in.fasta>.*
    inputBinding:
      position: 2
      prefix: "-6"

outputs:
  - id: output
    type: { type: array, items: File }
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

description: |
  Usage:   bwa index [options] <in.fasta>

  Options: -a STR    BWT construction algorithm: bwtsw or is [auto]
           -p STR    prefix of the index [same as fasta name]
           -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
           -6        index files named as <in.fasta>.64.* instead of <in.fasta>.*

  Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
           `-a div' do not work not for long genomes.


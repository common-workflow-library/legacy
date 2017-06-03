#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]
#TODO: Enable after this issue is fixed: https://github.com/common-workflow-language/cwltool/issues/80
#hints:
#  - $import: bwa-docker.yml

inputs:
  algorithm:
    type: string?
    inputBinding:
      prefix: -a
    doc: |
      BWT construction algorithm: bwtsw or is (Default: auto)
  sequences:
    type: File
    inputBinding:
      valueFrom: $(self.basename)
      position: 4
  block_size:
    type: int?
    inputBinding:
      prefix: -b
    doc: |
      Block size for the bwtsw algorithm (effective with -a bwtsw) (Default: 10000000)
baseCommand:
- bwa
- index

outputs:
  output:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
    outputBinding:
      glob: $(inputs.sequences.basename)


doc: |
  Usage:   bwa index [options] <in.fasta>

  Options: -a STR    BWT construction algorithm: bwtsw or is [auto]
           -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
           -6        index files named as <in.fasta>.64.* instead of <in.fasta>.*

  Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
           `-a div' do not work not for long genomes.


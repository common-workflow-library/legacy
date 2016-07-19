#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - $import: envvar-global.yml
  - class: InitialWorkDirRequirement
    listing:
     - $(inputs.reference_genome)

hints:
  - $import: bwa-docker.yml

inputs:
  reference_genome:
    type: File
    format: http://edamontology.org/format_1929  # FASTA
    inputBinding:
      position: 3

  algorithm:
    type: string?
    doc: |
      BWT construction algorithm: bwtsw or is (Default: auto)
    inputBinding:
      prefix: "-a"

  blocksize:
    type: int?
    doc: |
      Block size for the bwtsw algorithm (effective with -a bwtsw) (Default: 10000000)
    inputBinding:
      prefix: "-b"

outputs:
  output:
    type: File
    format: http://edamontology.org/format_1929  # FASTA
    secondaryFiles:
     - '.bai'
     - '.amb'
     - '.ann'
     - '.bwt'
     - '.pac'
     - '.sa'
    outputBinding:
      glob: $(inputs.reference_genome.basename)

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


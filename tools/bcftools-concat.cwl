#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

cwlVersion: "cwl:draft-3"

class: CommandLineTool

requirements:
  - $import: envvar-global.yml
  - $import: bcftools-docker.yml
  - class: InlineJavascriptRequirement

inputs:

  - id: filename
    type: string
    description: |
      Write output to a file [standard output]
    inputBinding:
      position: 1
      prefix: '-o'

  - id: vcfs
    type:
      type: array
      items: File
      secondaryFiles:
        - ".tbi"
    inputBinding:
      position: 2

  - id: allow_overlaps
    type:
      - 'null'
      - boolean
    description: >
      First coordinate of the next file can precede last record of the current file.
    inputBinding:
      position: 1
      prefix: '-a'

  - id: compact_PS
    type:
      - 'null'
      - boolean
    description: >
      Do not output PS tag at each site, only at the start of a new phase set block.
    inputBinding:
      position: 1
      prefix: '-c'

  - id: rm_dups
    type:
      - 'null'
      - string
    description: >
      Output duplicate records present in multiple files only once: <snps|indels|both|all|none>
    inputBinding:
      position: 1
      prefix: '-d'

  - id: remove_duplicates
    type:
      - 'null'
      - boolean
    description: |
      Alias for -d/--rm-dups none
    inputBinding:
      position: 1
      prefix: '-D'

  - id: file_list
    type:
      - 'null'
      - File
    description: |
      Read the list of files from a file.
    inputBinding:
      position: 1
      prefix: '-f'

  - id: ligate
    type:
      - 'null'
      - boolean
    description: >
      Ligate phased VCFs by matching phase at overlapping haplotypes
    inputBinding:
      position: 1
      prefix: '-l'

  - id: output_type
    type:
      - 'null'
      - string
    description: >
      <b|u|z|v> b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
    inputBinding:
      position: 1
      separate: false
      prefix: '-O'

  - id: min_PQ
    type:
      - 'null'
      - int
    description: >
      Break phase set if phasing quality is lower than <int> [30]
    inputBinding:
      position: 1
      prefix: '-q'

  - id: regions
    type:
      - 'null'
      - string
    description: |
      Restrict to comma-separated list of regions
    inputBinding:
      position: 1
      prefix: '-r'

  - id: regions_file
    type:
      - 'null'
      - File
    description: >
      Restrict to regions listed in a file
    inputBinding:
      position: 1
      prefix: '-R'

  - id: threads
    type:
      - 'null'
      - int
    description: >
      Number of extra output compression threads [0]
    inputBinding:
      position: 1
      prefix: '--threads'

outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.filename)

baseCommand:
  - bcftools
  - concat

description: |
  About:   Concatenate or combine VCF/BCF files. All source files must have the same sample
           columns appearing in the same order. The program can be used, for example, to
           concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel
           VCF into one. The input files must be sorted by chr and position. The files
           must be given in the correct order to produce sorted VCF on output unless
           the -a, --allow-overlaps option is specified.

  Usage:   bcftools concat [options] <A.vcf.gz> [<B.vcf.gz> [...]]


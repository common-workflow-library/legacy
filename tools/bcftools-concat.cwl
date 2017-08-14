#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: bcftools-docker.yml
- class: InlineJavascriptRequirement

inputs:
  allow_overlaps:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -a
    doc: 'First coordinate of the next file can precede last record of the current
      file.

      '
  remove_duplicates:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -D
    doc: |
      Alias for -d/--rm-dups none
  rm_dups:
    type: string?
    inputBinding:
      position: 1
      prefix: -d
    doc: 'Output duplicate records present in multiple files only once: <snps|indels|both|all|none>

      '
  ligate:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -l
    doc: 'Ligate phased VCFs by matching phase at overlapping haplotypes

      '
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: --threads
    doc: 'Number of extra output compression threads [0]

      '
  min_PQ:
    type: int?
    inputBinding:
      position: 1
      prefix: -q
    doc: 'Break phase set if phasing quality is lower than <int> [30]

      '
  filename:
    type: string
    inputBinding:
      position: 1
      prefix: -o
    doc: |
      Write output to a file [standard output]
  regions:
    type: string?
    inputBinding:
      position: 1
      prefix: -r
    doc: |
      Restrict to comma-separated list of regions
  file_list:
    type: File?
    inputBinding:
      position: 1
      prefix: -f
    doc: |
      Read the list of files from a file.
  compact_PS:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -c
    doc: 'Do not output PS tag at each site, only at the start of a new phase set
      block.

      '
  regions_file:
    type: File?
    inputBinding:
      position: 1
      prefix: -R
    doc: 'Restrict to regions listed in a file

      '
  vcfs:
    type:
      type: array
      items: File
    secondaryFiles:
      - .tbi
    inputBinding:
      position: 2

  output_type:
    type: string?
    inputBinding:
      position: 1
      separate: false
      prefix: -O
    doc: '<b|u|z|v> b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v:
      uncompressed VCF [v]

      '
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.filename)

baseCommand:
- bcftools
- concat

doc: |
  About:   Concatenate or combine VCF/BCF files. All source files must have the same sample
           columns appearing in the same order. The program can be used, for example, to
           concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel
           VCF into one. The input files must be sorted by chr and position. The files
           must be given in the correct order to produce sorted VCF on output unless
           the -a, --allow-overlaps option is specified.

  Usage:   bcftools concat [options] <A.vcf.gz> [<B.vcf.gz> [...]]


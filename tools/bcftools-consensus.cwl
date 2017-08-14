#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Dr. Barski Lab / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: bcftools-docker.yml
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement

inputs:
  vcf:
    type: File
    inputBinding:
      position: 2
    secondaryFiles:
    - .tbi

  reference:
    type: File
    inputBinding:
      position: 1
      prefix: -f
    doc: reference sequence in fasta format
  haplotype:
    type: int?
    inputBinding:
      position: 1
      prefix: -H
    doc: apply variants for the given haplotype <1|2>
  mask:
    type: File?
    inputBinding:
      position: 1
      prefix: -m
    doc: replace regions with N
  filename:
    type: string
    inputBinding:
      position: 1
      prefix: -o
    doc: write output to a file
  sample:
    type: string?
    inputBinding:
      position: 1
      prefix: -s
    doc: apply variants of the given sample
  iupac_codes:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -i
    doc: output variants in the form of IUPAC ambiguity codes
  chain:
    type: string?
    inputBinding:
      position: 1
      prefix: -c
    doc: write a chain file for liftover
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.filename)

  liftover:
    type: File
    outputBinding:
      glob: $(inputs.chain)

baseCommand: bcftools, consensus

stderr: ignore

doc: |
  About:   Create consensus sequence by applying VCF variants to a reference
           fasta file.

  Usage:   bcftools consensus [OPTIONS] <file.vcf>

  Options:
      -f, --fasta-ref <file>     reference sequence in fasta format
      -H, --haplotype <1|2>      apply variants for the given haplotype
      -i, --iupac-codes          output variants in the form of IUPAC ambiguity
  codes
      -m, --mask <file>          replace regions with N
      -o, --output <file>        write output to a file [standard output]
      -c, --chain <file>         write a chain file for liftover
      -s, --sample <name>        apply variants of the given sample

  Examples:

     # Get the consensus for one region. The fasta header lines are then expected
     # in the form ">chr:from-to".
     samtools faidx ref.fa 8:11870-11890 | bcftools consensus in.vcf.gz > out.fa


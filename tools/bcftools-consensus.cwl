#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Dr. Barski Lab / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

cwlVersion: 'cwl:draft-3'

class: CommandLineTool

requirements:
  - $import: envvar-global.yml
  - $import: bcftools-docker.yml
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement

inputs:

  - id: filename
    type: string
    description: |
      write output to a file
    inputBinding:
      position: 1
      prefix: '-o'

  - id: vcf
    type: File
    inputBinding:
      position: 2
    secondaryFiles:
      - .tbi

  - id: reference
    type: File
    description: |
      reference sequence in fasta format
    inputBinding:
      position: 1
      prefix: '-f'

  - id: haplotype
    type:
      - 'null'
      - int
    description: |
      apply variants for the given haplotype <1|2>
    inputBinding:
      position: 1
      prefix: '-H'

  - id: iupac_codes
    type:
      - 'null'
      - boolean
    description: >
      output variants in the form of IUPAC ambiguity codes
    inputBinding:
      position: 1
      prefix: '-i'

  - id: mask
    type:
      - 'null'
      - File
    description: |
      replace regions with N
    inputBinding:
      position: 1
      prefix: '-m'

  - id: chain
    type:
      - 'null'
      - string
    description: |
      write a chain file for liftover
    inputBinding:
      position: 1
      prefix: '-c'

  - id: sample
    type:
      - 'null'
      - string
    description: |
      apply variants of the given sample
    inputBinding:
      position: 1
      prefix: '-s'

outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.filename)

  - id: liftover
    type: File
    outputBinding:
      glob: $(inputs.chain)

baseCommand:
  - bcftools
  - consensus

arguments:
  - valueFrom: "2>/dev/null"
    position: 99
    shellQuote: false

description: |
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


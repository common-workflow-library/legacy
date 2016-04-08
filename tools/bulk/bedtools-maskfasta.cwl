#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: bedtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#fi'
    type:
      - 'null'
      - boolean
    description: |
      Input FASTA file
    inputBinding:
      position: 1
      prefix: '-fi'
  - id: '#bed'
    type:
      - 'null'
      - boolean
    description: |
      BED/GFF/VCF file of ranges to mask in -fi
    inputBinding:
      position: 1
      prefix: '-bed'
  - id: '#fo'
    type:
      - 'null'
      - boolean
    description: |
      Output FASTA file
    inputBinding:
      position: 1
      prefix: '-fo'
  - id: '#soft'
    type:
      - 'null'
      - boolean
    description: |
      Enforce "soft" masking.  That is, instead of masking with Ns,
      mask with lower-case bases.
    inputBinding:
      position: 1
      prefix: '-soft'
  - id: '#mc'
    type:
      - 'null'
      - boolean
    description: |
      Replace masking character.  That is, instead of masking
      with Ns, use another character.
    inputBinding:
      position: 1
      prefix: '-mc'
outputs:
  - id: '#stdoutfile'
    type: File
    outputBinding:
      glob:
        engine: 'cwl:JsonPointer'
        script: /job/stdoutfile
stdout:
  engine: 'cwl:JsonPointer'
  script: /job/stdoutfile
baseCommand:
  - bedtools
  - maskfasta
description: |
  Tool:    bedtools maskfasta (aka maskFastaFromBed)
  Version: v2.25.0
  Summary: Mask a fasta file based on feature coordinates.

  Usage:   bedtools maskfasta [OPTIONS] -fi <fasta> -fo <fasta> -bed <bed/gff/vcf>

  Options:
  	-fi	Input FASTA file
  	-bed	BED/GFF/VCF file of ranges to mask in -fi
  	-fo	Output FASTA file
  	-soft	Enforce "soft" masking.  That is, instead of masking with Ns,
  		mask with lower-case bases.
  	-mc	Replace masking character.  That is, instead of masking
  		with Ns, use another character.


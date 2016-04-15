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
      BED/GFF/VCF file of ranges to extract from -fi
    inputBinding:
      position: 1
      prefix: '-bed'
  - id: '#fo'
    type:
      - 'null'
      - boolean
    description: |
      Output file (can be FASTA or TAB-delimited)
    inputBinding:
      position: 1
      prefix: '-fo'
  - id: '#name'
    type:
      - 'null'
      - boolean
    description: |
      Use the name field for the FASTA header
    inputBinding:
      position: 1
      prefix: '-name'
  - id: '#split'
    type:
      - 'null'
      - boolean
    description: >
      given BED12 fmt., extract and concatenate the sequencesfrom the BED
      "blocks" (e.g., exons)
    inputBinding:
      position: 1
      prefix: '-split'
  - id: '#tab'
    type:
      - 'null'
      - boolean
    description: |
      Write output in TAB delimited format.
      - Default is FASTA format.
    inputBinding:
      position: 1
      prefix: '-tab'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Force strandedness. If the feature occupies the antisense,
      strand, the sequence will be reverse complemented.
      - By default, strand information is ignored.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#fullHeader'
    type:
      - 'null'
      - boolean
    description: |
      Use full fasta header.
      - By default, only the word before the first space or tab is used.
    inputBinding:
      position: 1
      prefix: '-fullHeader'
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
  - getfasta
description: |
  Tool:    bedtools getfasta (aka fastaFromBed)
  Version: v2.25.0
  Summary: Extract DNA sequences into a fasta file based on feature coordinates.

  Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta> 

  Options: 
  	-fi	Input FASTA file
  	-bed	BED/GFF/VCF file of ranges to extract from -fi
  	-fo	Output file (can be FASTA or TAB-delimited)
  	-name	Use the name field for the FASTA header
  	-split	given BED12 fmt., extract and concatenate the sequencesfrom the BED "blocks" (e.g., exons)
  	-tab	Write output in TAB delimited format.
  		- Default is FASTA format.

  	-s	Force strandedness. If the feature occupies the antisense,
  		strand, the sequence will be reverse complemented.
  		- By default, strand information is ignored.

  	-fullHeader	Use full fasta header.
  		- By default, only the word before the first space or tab is used.


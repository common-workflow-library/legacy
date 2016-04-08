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
  - id: '#i'
    type: File
    description: '<bam>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#bedpe'
    type:
      - 'null'
      - boolean
    description: |
      Write BEDPE format.
      - Requires BAM to be grouped or sorted by query.
    inputBinding:
      position: 1
      prefix: '-bedpe'
  - id: '#mate1'
    type:
      - 'null'
      - boolean
    description: "When writing BEDPE (-bedpe) format, \nalways report mate one as the first BEDPE \"block\".\n"
    inputBinding:
      position: 1
      prefix: '-mate1'
  - id: '#bed12'
    type:
      - 'null'
      - boolean
    description: |
      Write "blocked" BED format (aka "BED12"). Forces -split.
      http://genome-test.cse.ucsc.edu/FAQ/FAQformat#format1
    inputBinding:
      position: 1
      prefix: '-bed12'
  - id: '#split'
    type:
      - 'null'
      - boolean
    description: |
      Report "split" BAM alignments as separate BED entries.
      Splits only on N CIGAR operations.
    inputBinding:
      position: 1
      prefix: '-split'
  - id: '#splitD'
    type:
      - 'null'
      - boolean
    description: |
      Split alignments based on N and D CIGAR operators.
      Forces -split.
    inputBinding:
      position: 1
      prefix: '-splitD'
  - id: '#ed'
    type:
      - 'null'
      - boolean
    description: |
      Use BAM edit distance (NM tag) for BED score.
      - Default for BED is to use mapping quality.
      - Default for BEDPE is to use the minimum of
      the two mapping qualities for the pair.
      - When -ed is used with -bedpe, the total edit
      distance from the two mates is reported.
    inputBinding:
      position: 1
      prefix: '-ed'
  - id: '#tag'
    type:
      - 'null'
      - boolean
    description: |
      Use other NUMERIC BAM alignment tag for BED score.
      - Default for BED is to use mapping quality.
      Disallowed with BEDPE output.
    inputBinding:
      position: 1
      prefix: '-tag'
  - id: '#color'
    type:
      - 'null'
      - string
    description: |
      An R,G,B string for the color used with BED12 format.
      Default is (255,0,0).
    inputBinding:
      position: 1
      prefix: '-color'
  - id: '#cigar'
    type:
      - 'null'
      - string
    description: |
      Add the CIGAR string to the BED entry as a 7th column.
    inputBinding:
      position: 1
      prefix: '-cigar'
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
  - bamtobed
description: |
  Tool:    bedtools bamtobed (aka bamToBed)
  Version: v2.25.0
  Summary: Converts BAM alignments to BED6 or BEDPE format.

  Usage:   bedtools bamtobed [OPTIONS] -i <bam> 

  Options: 
  	-bedpe	Write BEDPE format.
  		- Requires BAM to be grouped or sorted by query.

  	-mate1	When writing BEDPE (-bedpe) format, 
  		always report mate one as the first BEDPE "block".

  	-bed12	Write "blocked" BED format (aka "BED12"). Forces -split.

  		http://genome-test.cse.ucsc.edu/FAQ/FAQformat#format1

  	-split	Report "split" BAM alignments as separate BED entries.
  		Splits only on N CIGAR operations.

  	-splitD	Split alignments based on N and D CIGAR operators.
  		Forces -split.

  	-ed	Use BAM edit distance (NM tag) for BED score.
  		- Default for BED is to use mapping quality.
  		- Default for BEDPE is to use the minimum of
  		  the two mapping qualities for the pair.
  		- When -ed is used with -bedpe, the total edit
  		  distance from the two mates is reported.

  	-tag	Use other NUMERIC BAM alignment tag for BED score.
  		- Default for BED is to use mapping quality.
  		  Disallowed with BEDPE output.

  	-color	An R,G,B string for the color used with BED12 format.
  		Default is (255,0,0).

  	-cigar	Add the CIGAR string to the BED entry as a 7th column.


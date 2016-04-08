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
  - id: '#g'
    type: boolean
    description: '<genome>'
    inputBinding:
      position: 3
      prefix: '-g'
  - id: '#i'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#mapq'
    type:
      - 'null'
      - int
    description: |
      Set the mappinq quality for the BAM records.
      (INT) Default: 255
    inputBinding:
      position: 1
      prefix: '-mapq'
  - id: '#bed12'
    type:
      - 'null'
      - string
    description: |
      The BED file is in BED12 format.  The BAM CIGAR
      string will reflect BED "blocks".
    inputBinding:
      position: 1
      prefix: '-bed12'
  - id: '#ubam'
    type:
      - 'null'
      - boolean
    description: |
      Write uncompressed BAM output. Default writes compressed BAM.
    inputBinding:
      position: 1
      prefix: '-ubam'
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
  - bedtobam
description: |
  Tool:    bedtools bedtobam (aka bedToBam)
  Version: v2.25.0
  Summary: Converts feature records to BAM format.

  Usage:   bedtools bedtobam [OPTIONS] -i <bed/gff/vcf> -g <genome>

  Options: 
  	-mapq	Set the mappinq quality for the BAM records.
  		(INT) Default: 255

  	-bed12	The BED file is in BED12 format.  The BAM CIGAR
  		string will reflect BED "blocks".

  	-ubam	Write uncompressed BAM output. Default writes compressed BAM.

  Notes: 
  	(1)  BED files must be at least BED4 to create BAM (needs name field).


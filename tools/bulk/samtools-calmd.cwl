#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: samtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#reffasta'
    type: File
    description: '<ref.fasta>'
    inputBinding:
      position: 4
  - id: '#alnbam'
    type: File
    description: '<aln.bam>'
    inputBinding:
      position: 3
  - id: '#e'
    type:
      - 'null'
      - boolean
    description: |
      change identical bases to '='
    inputBinding:
      position: 1
      prefix: '-e'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: "      uncompressed BAM output (for piping)\n"
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#b'
    type:
      - 'null'
      - boolean
    description: "      compressed BAM output\n"
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: "      the input is SAM with header\n"
    inputBinding:
      position: 1
      prefix: '-S'
  - id: '#A'
    type:
      - 'null'
      - string
    description: "      modify the quality string\n"
    inputBinding:
      position: 1
      prefix: '-A'
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: "      compute the BQ tag (without -A) or cap baseQ by BAQ (with -A)\n"
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#E'
    type:
      - 'null'
      - boolean
    description: "      extended BAQ for better sensitivity but lower specificity\n"
    inputBinding:
      position: 1
      prefix: '-E'
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
  - samtools
  - calmd
description: |-
  Usage:   samtools calmd [-eubrS] <aln.bam> <ref.fasta>

  Options: -e       change identical bases to '='
           -u       uncompressed BAM output (for piping)
           -b       compressed BAM output
           -S       the input is SAM with header
           -A       modify the quality string
           -r       compute the BQ tag (without -A) or cap baseQ by BAQ (with -A)
           -E       extended BAQ for better sensitivity but lower specificity


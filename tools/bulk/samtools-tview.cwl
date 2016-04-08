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
  - id: '#[reffasta]'
    type:
      - 'null'
      - File
    description: '[ref.fasta]'
    inputBinding:
      position: 3
  - id: '#alnbam'
    type: File
    description: '<aln.bam>'
    inputBinding:
      position: 2
  - id: '#d'
    type:
      - 'null'
      - boolean
    description: "display      output as (H)tml or (C)urses or (T)ext \n"
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#p'
    type:
      - 'null'
      - boolean
    description: |
      chr:pos      go directly to this position
    inputBinding:
      position: 1
      prefix: '-p'
  - id: '#s'
    type:
      - 'null'
      - string
    description: |
      STR          display only reads from this sample or group
    inputBinding:
      position: 1
      prefix: '-s'
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
  - tview
description: |
  Usage: samtools tview [options] <aln.bam> [ref.fasta]
  Options:
     -d display      output as (H)tml or (C)urses or (T)ext 
     -p chr:pos      go directly to this position
     -s STR          display only reads from this sample or group


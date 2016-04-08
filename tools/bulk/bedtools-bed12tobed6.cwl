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
    type: boolean
    description: '<bed12>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#n'
    type:
      - 'null'
      - boolean
    description: |
      Force the score to be the (1-based) block number from the BED12.
    inputBinding:
      position: 1
      prefix: '-n'
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
  - bed12tobed6
description: |
  Tool:    bedtools bed12tobed6 (aka bed12ToBed6)
  Version: v2.25.0
  Summary: Splits BED12 features into discrete BED6 features.

  Usage:   bedtools bed12tobed6 [OPTIONS] -i <bed12>

  Options: 
  	-n	Force the score to be the (1-based) block number from the BED12.


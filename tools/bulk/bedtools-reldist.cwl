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
  - id: '#b'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 3
      prefix: '-b'
  - id: '#a'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-a'
  - id: '#detail'
    type:
      - 'null'
      - boolean
    description: "Instead of a summary, report the relative\t\t distance for each interval in A\n"
    inputBinding:
      position: 1
      prefix: '-detail'
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
  - reldist
description: |
  Tool:    bedtools reldist
  Version: v2.25.0
  Summary: Calculate the relative distance distribution b/w two feature files.

  Usage:   bedtools reldist [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

  Options: 
  	-detail	Instead of a summary, report the relative		 distance for each interval in A


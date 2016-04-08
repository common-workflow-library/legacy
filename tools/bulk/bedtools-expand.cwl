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
    type:
      - 'null'
      - boolean
    description: |
      Input file. Assumes "stdin" if omitted.
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: "\tSpecify the column (1-based) that should be summarized.\n- Required.\n"
    inputBinding:
      position: 1
      prefix: '-c'
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
  - expand
description: |
  Tool:    bedtools expand 
  Version: v2.25.0
  Summary: Replicate lines in a file based on columns of comma-separated values.

  Usage:	 bedtools expand -c [COLS] 
  Options: 
  	-i	Input file. Assumes "stdin" if omitted.

  	-c 	Specify the column (1-based) that should be summarized.
  		- Required.
  Examples: 
    $ cat test.txt
    chr1	10	20	1,2,3	10,20,30
    chr1	40	50	4,5,6	40,50,60

    $ bedtools expand test.txt -c 5
    chr1	10	20	1,2,3	10
    chr1	10	20	1,2,3	20
    chr1	10	20	1,2,3	30
    chr1	40	50	4,5,6	40
    chr1	40	50	4,5,6	50
    chr1	40	50	4,5,6	60

    $ bedtools expand test.txt -c 4,5
    chr1	10	20	1	10
    chr1	10	20	2	20
    chr1	10	20	3	30
    chr1	40	50	4	40
    chr1	40	50	5	50
    chr1	40	50	6	60


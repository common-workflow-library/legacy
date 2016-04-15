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
      position: 2
      prefix: '-g'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      The length of the intervals to generate.
      - Default = 100.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#n'
    type:
      - 'null'
      - int
    description: |
      The number of intervals to generate.
      - Default = 1,000,000.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#seed'
    type:
      - 'null'
      - int
    description: |
      Supply an integer seed for the shuffling.
      - By default, the seed is chosen automatically.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-seed'
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
  - random
description: |
  Tool:    bedtools random (aka randomBed)
  Version: v2.25.0
  Summary: Generate random intervals among a genome.

  Usage:   bedtools random [OPTIONS] -g <genome>

  Options: 
  	-l	The length of the intervals to generate.
  		- Default = 100.
  		- (INTEGER)

  	-n	The number of intervals to generate.
  		- Default = 1,000,000.
  		- (INTEGER)

  	-seed	Supply an integer seed for the shuffling.
  		- By default, the seed is chosen automatically.
  		- (INTEGER)

  Notes: 
  	(1)  The genome file should tab delimited and structured as follows:
  	     <chromName><TAB><chromSize>

  	For example, Human (hg19):
  	chr1	249250621
  	chr2	243199373
  	...
  	chr18_gl000207_random	4262

  Tips: 
  	One can use the UCSC Genome Browser's MySQL database to extract
  	chromosome sizes. For example, H. sapiens:

  	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
  	"select chrom, size from hg19.chromInfo"  > hg19.genome


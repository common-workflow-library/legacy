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
    description: '<BEDPE>'
    inputBinding:
      position: 3
      prefix: '-b'
  - id: '#a'
    type: File
    description: '<BEDPE>'
    inputBinding:
      position: 2
      prefix: '-a'
  - id: '#f'
    type:
      - 'null'
      - boolean
    description: |
      Minimum overlap required as fraction of A (e.g. 0.05).
      Default is 1E-9 (effectively 1bp).
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#type'
    type:
      - 'null'
      - boolean
    description: |
      	Approach to reporting overlaps between A and B.
      neither	Report overlaps if neither end of A overlaps B.
      either	Report overlaps if either ends of A overlap B.
      both	Report overlaps if both ends of A overlap B.
      notboth	Report overlaps if one or neither of A's overlap B.
      - Default = both.
    inputBinding:
      position: 1
      prefix: '-type'
  - id: '#slop'
    type:
      - 'null'
      - boolean
    description: |
      	The amount of slop (in b.p.). to be added to each footprint of A.
      *Note*: Slop is subtracted from start1 and start2
      and added to end1 and end2.
      - Default = 0.
    inputBinding:
      position: 1
      prefix: '-slop'
  - id: '#ss'
    type:
      - 'null'
      - boolean
    description: |
      Add slop based to each BEDPE footprint based on strand.
      - If strand is "+", slop is only added to the end coordinates.
      - If strand is "-", slop is only added to the start coordinates.
      - By default, slop is added in both directions.
    inputBinding:
      position: 1
      prefix: '-ss'
  - id: '#is'
    type:
      - 'null'
      - boolean
    description: |
      Ignore strands when searching for overlaps.
      - By default, strands are enforced.
    inputBinding:
      position: 1
      prefix: '-is'
  - id: '#rdn'
    type:
      - 'null'
      - boolean
    description: |
      Require the hits to have different names (i.e. avoid self-hits).
      - By default, same names are allowed.
      Refer to the BEDTools manual for BEDPE format.
    inputBinding:
      position: 1
      prefix: '-rdn'
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
  - pairtopair
description: |
  Tool:    bedtools pairtopair (aka pairToPair)
  Version: v2.25.0
  Summary: Report overlaps between two paired-end BED files (BEDPE).

  Usage:   bedtools pairtopair [OPTIONS] -a <BEDPE> -b <BEDPE>

  Options: 
  	-f	Minimum overlap required as fraction of A (e.g. 0.05).
  		Default is 1E-9 (effectively 1bp).

  	-type 	Approach to reporting overlaps between A and B.

  		neither	Report overlaps if neither end of A overlaps B.
  		either	Report overlaps if either ends of A overlap B.
  		both	Report overlaps if both ends of A overlap B.
  		notboth	Report overlaps if one or neither of A's overlap B.
  		- Default = both.

  	-slop 	The amount of slop (in b.p.). to be added to each footprint of A.
  		*Note*: Slop is subtracted from start1 and start2
  			and added to end1 and end2.

  		- Default = 0.

  	-ss	Add slop based to each BEDPE footprint based on strand.
  		- If strand is "+", slop is only added to the end coordinates.
  		- If strand is "-", slop is only added to the start coordinates.
  		- By default, slop is added in both directions.

  	-is	Ignore strands when searching for overlaps.
  		- By default, strands are enforced.

  	-rdn	Require the hits to have different names (i.e. avoid self-hits).
  		- By default, same names are allowed.

  Refer to the BEDTools manual for BEDPE format.


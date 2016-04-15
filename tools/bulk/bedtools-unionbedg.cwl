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
  - id: '#FILEn'
    type: boolean
    description: FILEn
    inputBinding:
      position: 5
  - id: '#'
    type: boolean
    description: ..
    inputBinding:
      position: 4
  - id: '#FILE2'
    type: File
    description: FILE2
    inputBinding:
      position: 3
  - id: '#i'
    type: File
    description: FILE1
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: "\tPrint a header line.\n(chrom/start/end + names of each file).\n"
    inputBinding:
      position: 1
      prefix: '-header'
  - id: '#names'
    type:
      - 'null'
      - boolean
    description: "\tA list of names (one/file) to describe each file in -i.\nThese names will be printed in the header line.\n"
    inputBinding:
      position: 1
      prefix: '-names'
  - id: '#g'
    type:
      - 'null'
      - string
    description: "\tUse genome file to calculate empty regions.\n- STRING.\n"
    inputBinding:
      position: 1
      prefix: '-g'
  - id: '#empty'
    type:
      - 'null'
      - File
    description: |
      	Report empty regions (i.e., start/end intervals w/o
      values in all files).
      - Requires the '-g FILE' parameter.
    inputBinding:
      position: 1
      prefix: '-empty'
  - id: '#filler'
    type:
      - 'null'
      - string
    description: "TEXT\tUse TEXT when representing intervals having no value.\n- Default is '0', but you can use 'N/A' or any text.\n"
    inputBinding:
      position: 1
      prefix: '-filler'
  - id: '#examples'
    type:
      - 'null'
      - boolean
    description: |
      Show detailed usage examples.
    inputBinding:
      position: 1
      prefix: '-examples'
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
  - unionbedg
description: |
  Tool:    bedtools unionbedg (aka unionBedGraphs)
  Version: v2.25.0
  Summary: Combines multiple BedGraph files into a single file,
  	 allowing coverage comparisons between them.

  Usage:   bedtools unionbedg [OPTIONS] -i FILE1 FILE2 .. FILEn
  	 Assumes that each BedGraph file is sorted by chrom/start 
  	 and that the intervals in each are non-overlapping.

  Options: 
  	-header		Print a header line.
  			(chrom/start/end + names of each file).

  	-names		A list of names (one/file) to describe each file in -i.
  			These names will be printed in the header line.

  	-g		Use genome file to calculate empty regions.
  			- STRING.

  	-empty		Report empty regions (i.e., start/end intervals w/o
  			values in all files).
  			- Requires the '-g FILE' parameter.

  	-filler TEXT	Use TEXT when representing intervals having no value.
  			- Default is '0', but you can use 'N/A' or any text.

  	-examples	Show detailed usage examples.


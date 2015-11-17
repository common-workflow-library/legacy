#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: bedtools-docker.cwl
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#cluster'
    type:
      - 'null'
      - boolean
    description: |
      Invoke Ryan Layers's clustering algorithm.
    inputBinding:
      position: 1
      prefix: '-cluster'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: |
      Print a header line.
      (chrom/start/end + names of each file).
    inputBinding:
      position: 1
      prefix: '-header'
  - id: '#names'
    type:
      - 'null'
      - boolean
    description: |
      A list of names (one/file) to describe each file in -i.
      These names will be printed in the header line.
    inputBinding:
      position: 1
      prefix: '-names'
  - id: '#g'
    type:
      - 'null'
      - string
    description: |
      Use genome file to calculate empty regions.
      - STRING.
    inputBinding:
      position: 1
      prefix: '-g'
  - id: '#empty'
    type:
      - 'null'
      - boolean
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
  - multiinter
description: |

  Tool:    bedtools multiinter (aka multiIntersectBed)
  Version: v2.25.0
  Summary: Identifies common intervals among multiple
  	 BED/GFF/VCF files.

  Usage:   bedtools multiinter [OPTIONS] -i FILE1 FILE2 .. FILEn
  	 Requires that each interval file is sorted by chrom/start. 

  Options: 
  	-cluster	Invoke Ryan Layers's clustering algorithm.

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


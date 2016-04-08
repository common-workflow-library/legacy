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
  - id: '#s'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#e'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-e'
  - id: '#l'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#f'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#i'
    type:
      - 'null'
      - File
    description: "--input (file)\tBED input file (req'd).\n"
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#n'
    type:
      - 'null'
      - int
    description: "--number (int)\tNumber of files to create (req'd).\n"
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#p'
    type:
      - 'null'
      - string
    description: "--prefix (string)\tOutput BED file prefix.\n"
    inputBinding:
      position: 1
      prefix: '-p'
  - id: '#a'
    type:
      - 'null'
      - string
    description: |
      --algorithm (string) Algorithm used to split data.
      * size (default): uses a heuristic algorithm to group the items
      so all files contain the ~ same number of bases
      * simple : route records such that each split file has
      approximately equal records (like Unix split).
    inputBinding:
      position: 1
      prefix: '-a'
  - id: '#h'
    type:
      - 'null'
      - boolean
    description: "--help\t\tPrint help (this screen).\n"
    inputBinding:
      position: 1
      prefix: '-h'
  - id: '#v'
    type:
      - 'null'
      - boolean
    description: "--version\t\tPrint version.\nNote: This programs stores the input BED records in memory.\n"
    inputBinding:
      position: 1
      prefix: '-v'
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
  - split
description: |
  Error: num_chunks==0.


  Tool:    bedtools split
  Version: v2.25.0
  Summary: Split a Bed file.

  Usage:   bedtools split [OPTIONS] -i <bed> -n number-of-files

  Options: 
  	-i|--input (file)	BED input file (req'd).
  	-n|--number (int)	Number of files to create (req'd).
  	-p|--prefix (string)	Output BED file prefix.
  	-a|--algorithm (string) Algorithm used to split data.
  		* size (default): uses a heuristic algorithm to group the items
  		  so all files contain the ~ same number of bases
  		* simple : route records such that each split file has
  		  approximately equal records (like Unix split).

  	-h|--help		Print help (this screen).
  	-v|--version		Print version.


  Note: This programs stores the input BED records in memory.


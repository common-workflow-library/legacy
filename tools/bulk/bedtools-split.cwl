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
  - id: '#i'
    type: File
    description: |
      <bed>
    inputBinding:
      position: 2
      prefix: '-i'
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



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
      Input file. Use "stdin" for pipes.
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#cols'
    type:
      - 'null'
      - boolean
    description: |
      Specify the columns (1-based) for the starts and ends of the
      features for which you'd like to compute the overlap/distance.
      The columns must be listed in the following order:
      start1,end1,start2,end2
      Example:
      $ windowBed -a A.bed -b B.bed -w 10
      chr1 10  20  A   chr1    15  25  B
      chr1 10  20  C   chr1    25  35  D
      $ windowBed -a A.bed -b B.bed -w 10 | overlap -i stdin -cols 2,3,6,7
      chr1 10  20  A   chr1    15  25  B   5
      chr1 10  20  C   chr1    25  35  D   -5
    inputBinding:
      position: 1
      prefix: '-cols'
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
  - overlap
description: |
  Tool:    bedtools overlap (aka getOverlap)
  Version: v2.25.0
  Summary: Computes the amount of overlap (positive values)
  	 or distance (negative values) between genome features
  	 and reports the result at the end of the same line.

  Options: 
  	-i	Input file. Use "stdin" for pipes.

  	-cols	Specify the columns (1-based) for the starts and ends of the
  		features for which you'd like to compute the overlap/distance.
  		The columns must be listed in the following order: 

  		start1,end1,start2,end2

  Example: 
  	$ windowBed -a A.bed -b B.bed -w 10
  	chr1 10  20  A   chr1    15  25  B
  	chr1 10  20  C   chr1    25  35  D

  	$ windowBed -a A.bed -b B.bed -w 10 | overlap -i stdin -cols 2,3,6,7
  	chr1 10  20  A   chr1    15  25  B   5
  	chr1 10  20  C   chr1    25  35  D   -5


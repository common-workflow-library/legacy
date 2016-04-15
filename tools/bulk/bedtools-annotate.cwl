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
  - id: '#FILE2FILEn'
    type: File
    description: FILE2..FILEn
    inputBinding:
      position: 4
  - id: '#files'
    type: File
    description: FILE1
    inputBinding:
      position: 3
      prefix: '-files'
  - id: '#i'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#names'
    type:
      - 'null'
      - boolean
    description: |
      A list of names (one / file) to describe each file in -i.
      These names will be printed as a header line.
    inputBinding:
      position: 1
      prefix: '-names'
  - id: '#counts'
    type:
      - 'null'
      - boolean
    description: |
      Report the count of features in each file that overlap -i.
      - Default is to report the fraction of -i covered by each file.
    inputBinding:
      position: 1
      prefix: '-counts'
  - id: '#both'
    type:
      - 'null'
      - boolean
    description: |
      Report the counts followed by the % coverage.
      - Default is to report the fraction of -i covered by each file.
    inputBinding:
      position: 1
      prefix: '-both'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Require same strandedness.  That is, only counts overlaps
      on the _same_ strand.
      - By default, overlaps are counted without respect to strand.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: |
      Require different strandedness.  That is, only count overlaps
      on the _opposite_ strand.
      - By default, overlaps are counted without respect to strand.
    inputBinding:
      position: 1
      prefix: '-S'
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
  - annotate
description: |
  Tool:    bedtools annotate (aka annotateBed)
  Version: v2.25.0
  Summary: Annotates the depth & breadth of coverage of features from mult. files
  	 on the intervals in -i.

  Usage:   bedtools annotate [OPTIONS] -i <bed/gff/vcf> -files FILE1 FILE2..FILEn

  Options: 
  	-names	A list of names (one / file) to describe each file in -i.
  		These names will be printed as a header line.

  	-counts	Report the count of features in each file that overlap -i.
  		- Default is to report the fraction of -i covered by each file.

  	-both	Report the counts followed by the % coverage.
  		- Default is to report the fraction of -i covered by each file.

  	-s	Require same strandedness.  That is, only counts overlaps
  		on the _same_ strand.
  		- By default, overlaps are counted without respect to strand.

  	-S	Require different strandedness.  That is, only count overlaps
  		on the _opposite_ strand.
  		- By default, overlaps are counted without respect to strand.


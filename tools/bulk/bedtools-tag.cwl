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
  - id: '#LABn'
    type: boolean
    description: LABn
    inputBinding:
      position: 8
  - id: '#labels'
    type: boolean
    description: LAB1
    inputBinding:
      position: 6
      prefix: '-labels'
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
  - id: '#files'
    type: File
    description: FILE1
    inputBinding:
      position: 3
      prefix: '-files'
  - id: '#i'
    type: File
    description: '<BAM>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: >
      Require overlaps on the same strand.  That is, only tag alignments that
      have the same

      strand as a feature in the annotation file(s).
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: >
      Require overlaps on the opposite strand.  That is, only tag alignments that
      have the opposite

      strand as a feature in the annotation file(s).
    inputBinding:
      position: 1
      prefix: '-S'
  - id: '#f'
    type:
      - 'null'
      - float
    description: |
      Minimum overlap required as a fraction of the alignment.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#tag'
    type:
      - 'null'
      - string
    description: |
      Dictate what the tag should be. Default is YB.
      - STRING (two characters, e.g., YK)
    inputBinding:
      position: 1
      prefix: '-tag'
  - id: '#names'
    type:
      - 'null'
      - boolean
    description: |
      Use the name field from the annotation files to populate tags.
      By default, the -labels values are used.
    inputBinding:
      position: 1
      prefix: '-names'
  - id: '#scores'
    type:
      - 'null'
      - boolean
    description: |
      Use the score field from the annotation files to populate tags.
      By default, the -labels values are used.
    inputBinding:
      position: 1
      prefix: '-scores'
  - id: '#intervals'
    type:
      - 'null'
      - boolean
    description: >
      Use the full interval (including name, score, and strand) to populate
      tags.

      Requires the -labels option to identify from which file the interval came.
    inputBinding:
      position: 1
      prefix: '-intervals'
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
  - tag
description: |
  Tool:    bedtools tag (aka tagBam)
  Version: v2.25.0
  Summary: Annotates a BAM file based on overlaps with multiple BED/GFF/VCF files
  	 on the intervals in -i.

  Usage:   bedtools tag [OPTIONS] -i <BAM> -files FILE1 .. FILEn  -labels LAB1 .. LABn

  Options: 
  	-s	Require overlaps on the same strand.  That is, only tag alignments that have the same
  		strand as a feature in the annotation file(s).

  	-S	Require overlaps on the opposite strand.  That is, only tag alignments that have the opposite
  		strand as a feature in the annotation file(s).

  	-f	Minimum overlap required as a fraction of the alignment.
  		- Default is 1E-9 (i.e., 1bp).
  		- FLOAT (e.g. 0.50)

  	-tag	Dictate what the tag should be. Default is YB.
  		- STRING (two characters, e.g., YK)

  	-names	Use the name field from the annotation files to populate tags.
  		By default, the -labels values are used.

  	-scores	Use the score field from the annotation files to populate tags.
  		By default, the -labels values are used.

  	-intervals	Use the full interval (including name, score, and strand) to populate tags.
  			Requires the -labels option to identify from which file the interval came.


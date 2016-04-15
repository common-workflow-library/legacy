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
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#path'
    type:
      - 'null'
      - string
    description: |
      The full path to which the IGV snapshots should be written.
      (STRING) Default: ./
    inputBinding:
      position: 1
      prefix: '-path'
  - id: '#sess'
    type:
      - 'null'
      - string
    description: |
      The full path to an existing IGV session file to be 
      loaded prior to taking snapshots.
      (STRING) Default is for no session to be loaded.
    inputBinding:
      position: 1
      prefix: '-sess'
  - id: '#sort'
    type:
      - 'null'
      - boolean
    description: |
      The type of BAM sorting you would like to apply to each image. 
      Options: base, position, strand, quality, sample, and readGroup
      Default is to apply no sorting at all.
    inputBinding:
      position: 1
      prefix: '-sort'
  - id: '#clps'
    type:
      - 'null'
      - boolean
    description: "Collapse the aligned reads prior to taking a snapshot. \nDefault is to no collapse.\n"
    inputBinding:
      position: 1
      prefix: '-clps'
  - id: '#name'
    type:
      - 'null'
      - boolean
    description: "Use the \"name\" field (column 4) for each image's filename. \nDefault is to use the \"chr:start-pos.ext\".\n"
    inputBinding:
      position: 1
      prefix: '-name'
  - id: '#slop'
    type:
      - 'null'
      - int
    description: |
      Number of flanking base pairs on the left & right of the image.
      - (INT) Default = 0.
    inputBinding:
      position: 1
      prefix: '-slop'
  - id: '#img'
    type:
      - 'null'
      - boolean
    description: "The type of image to be created. \nOptions: png, eps, svg\nDefault is png.\n"
    inputBinding:
      position: 1
      prefix: '-img'
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
  - igv
description: |
  Tool:    bedtools igv (aka bedToIgv)
  Version: v2.25.0
  Summary: Creates a batch script to create IGV images 
           at each interval defined in a BED/GFF/VCF file.

  Usage:   bedtools igv [OPTIONS] -i <bed/gff/vcf>

  Options: 
  	-path	The full path to which the IGV snapshots should be written.
  		(STRING) Default: ./

  	-sess	The full path to an existing IGV session file to be 
  		loaded prior to taking snapshots.

  		(STRING) Default is for no session to be loaded.

  	-sort	The type of BAM sorting you would like to apply to each image. 
  		Options: base, position, strand, quality, sample, and readGroup
  		Default is to apply no sorting at all.

  	-clps	Collapse the aligned reads prior to taking a snapshot. 
  		Default is to no collapse.

  	-name	Use the "name" field (column 4) for each image's filename. 
  		Default is to use the "chr:start-pos.ext".

  	-slop	Number of flanking base pairs on the left & right of the image.
  		- (INT) Default = 0.

  	-img	The type of image to be created. 
  		Options: png, eps, svg
  		Default is png.

  Notes: 
  	(1)  The resulting script is meant to be run from within IGV.
  	(2)  Unless you use the -sess option, it is assumed that prior to 
  		running the script, you've loaded the proper genome and tracks.


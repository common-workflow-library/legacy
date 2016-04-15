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
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 3
      prefix: '-b'
  - id: '#a'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-a'
  - id: '#abam'
    type:
      - 'null'
      - boolean
    description: >
      The A input file is in BAM format.  Output will be BAM as well. Replaces
      -a.
    inputBinding:
      position: 1
      prefix: '-abam'
  - id: '#ubam'
    type:
      - 'null'
      - boolean
    description: |
      Write uncompressed BAM output. Default writes compressed BAM.
    inputBinding:
      position: 1
      prefix: '-ubam'
  - id: '#bed'
    type:
      - 'null'
      - boolean
    description: |
      When using BAM input (-abam), write output as BED. The default
      is to write output in BAM when using -abam.
    inputBinding:
      position: 1
      prefix: '-bed'
  - id: '#w'
    type:
      - 'null'
      - int
    description: |
      Base pairs added upstream and downstream of each entry
      in A when searching for overlaps in B.
      - Creates symterical "windows" around A.
      - Default is 1000 bp.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-w'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      Base pairs added upstream (left of) of each entry
      in A when searching for overlaps in B.
      - Allows one to define assymterical "windows".
      - Default is 1000 bp.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#r'
    type:
      - 'null'
      - int
    description: |
      Base pairs added downstream (right of) of each entry
      in A when searching for overlaps in B.
      - Allows one to define assymterical "windows".
      - Default is 1000 bp.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#sw'
    type:
      - 'null'
      - boolean
    description: |
      Define -l and -r based on strand.  For example if used, -l 500
      for a negative-stranded feature will add 500 bp downstream.
      - Default = disabled.
    inputBinding:
      position: 1
      prefix: '-sw'
  - id: '#sm'
    type:
      - 'null'
      - boolean
    description: |
      Only report hits in B that overlap A on the _same_ strand.
      - By default, overlaps are reported without respect to strand.
    inputBinding:
      position: 1
      prefix: '-sm'
  - id: '#Sm'
    type:
      - 'null'
      - boolean
    description: |
      Only report hits in B that overlap A on the _opposite_ strand.
      - By default, overlaps are reported without respect to strand.
    inputBinding:
      position: 1
      prefix: '-Sm'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: |
      Write the original A entry _once_ if _any_ overlaps found in B.
      - In other words, just report the fact >=1 hit was found.
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: |
      For each entry in A, report the number of overlaps with B.
      - Reports 0 for A entries that have no overlap with B.
      - Overlaps restricted by -f.
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#v'
    type:
      - 'null'
      - boolean
    description: |
      Only report those entries in A that have _no overlaps_ with B.
      - Similar to "grep -v."
    inputBinding:
      position: 1
      prefix: '-v'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: |
      Print the header from the A file prior to results.
    inputBinding:
      position: 1
      prefix: '-header'
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
  - window
description: |
  Tool:    bedtools window (aka windowBed)
  Version: v2.25.0
  Summary: Examines a "window" around each feature in A and
  	 reports all features in B that overlap the window. For each
  	 overlap the entire entry in A and B are reported.

  Usage:   bedtools window [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

  Options: 
  	-abam	The A input file is in BAM format.  Output will be BAM as well. Replaces -a.

  	-ubam	Write uncompressed BAM output. Default writes compressed BAM.

  	-bed	When using BAM input (-abam), write output as BED. The default
  		is to write output in BAM when using -abam.

  	-w	Base pairs added upstream and downstream of each entry
  		in A when searching for overlaps in B.
  		- Creates symterical "windows" around A.
  		- Default is 1000 bp.
  		- (INTEGER)

  	-l	Base pairs added upstream (left of) of each entry
  		in A when searching for overlaps in B.
  		- Allows one to define assymterical "windows".
  		- Default is 1000 bp.
  		- (INTEGER)

  	-r	Base pairs added downstream (right of) of each entry
  		in A when searching for overlaps in B.
  		- Allows one to define assymterical "windows".
  		- Default is 1000 bp.
  		- (INTEGER)

  	-sw	Define -l and -r based on strand.  For example if used, -l 500
  		for a negative-stranded feature will add 500 bp downstream.
  		- Default = disabled.

  	-sm	Only report hits in B that overlap A on the _same_ strand.
  		- By default, overlaps are reported without respect to strand.

  	-Sm	Only report hits in B that overlap A on the _opposite_ strand.
  		- By default, overlaps are reported without respect to strand.

  	-u	Write the original A entry _once_ if _any_ overlaps found in B.
  		- In other words, just report the fact >=1 hit was found.

  	-c	For each entry in A, report the number of overlaps with B.
  		- Reports 0 for A entries that have no overlap with B.
  		- Overlaps restricted by -f.

  	-v	Only report those entries in A that have _no overlaps_ with B.
  		- Similar to "grep -v."

  	-header	Print the header from the A file prior to results.


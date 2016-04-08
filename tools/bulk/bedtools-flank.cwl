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
  - id: '#g'
    type: boolean
    description: '<genome>'
    inputBinding:
      position: 3
      prefix: '-g'
  - id: '#i'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#b'
    type:
      - 'null'
      - float
    description: |
      Create flanking interval(s) using -b base pairs in each direction.
      - (Integer) or (Float, e.g. 0.1) if used with -pct.
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#l'
    type:
      - 'null'
      - float
    description: |
      The number of base pairs that a flank should start from
      orig. start coordinate.
      - (Integer) or (Float, e.g. 0.1) if used with -pct.
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#r'
    type:
      - 'null'
      - float
    description: |
      The number of base pairs that a flank should end from
      orig. end coordinate.
      - (Integer) or (Float, e.g. 0.1) if used with -pct.
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Define -l and -r based on strand.
      E.g. if used, -l 500 for a negative-stranded feature,
      it will start the flank 500 bp downstream.  Default = false.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#pct'
    type:
      - 'null'
      - boolean
    description: |
      Define -l and -r as a fraction of the feature's length.
      E.g. if used on a 1000bp feature, -l 0.50,
      will add 500 bp "upstream".  Default = false.
    inputBinding:
      position: 1
      prefix: '-pct'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: |
      Print the header from the input file prior to results.
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
  - flank
description: |
  Tool:    bedtools flank (aka flankBed)
  Version: v2.25.0
  Summary: Creates flanking interval(s) for each BED/GFF/VCF feature.

  Usage:   bedtools flank [OPTIONS] -i <bed/gff/vcf> -g <genome> [-b <int> or (-l and -r)]

  Options: 
  	-b	Create flanking interval(s) using -b base pairs in each direction.
  		- (Integer) or (Float, e.g. 0.1) if used with -pct.

  	-l	The number of base pairs that a flank should start from
  		orig. start coordinate.
  		- (Integer) or (Float, e.g. 0.1) if used with -pct.

  	-r	The number of base pairs that a flank should end from
  		orig. end coordinate.
  		- (Integer) or (Float, e.g. 0.1) if used with -pct.

  	-s	Define -l and -r based on strand.
  		E.g. if used, -l 500 for a negative-stranded feature, 
  		it will start the flank 500 bp downstream.  Default = false.

  	-pct	Define -l and -r as a fraction of the feature's length.
  		E.g. if used on a 1000bp feature, -l 0.50, 
  		will add 500 bp "upstream".  Default = false.

  	-header	Print the header from the input file prior to results.

  Notes: 
  	(1)  Starts will be set to 0 if options would force it below 0.
  	(2)  Ends will be set to the chromosome length if requested flank would
  	force it above the max chrom length.
  	(3)  In contrast to slop, which _extends_ intervals, bedtools flank
  	creates new intervals from the regions just up- and down-stream
  	of your existing intervals.
  	(4)  The genome file should tab delimited and structured as follows:

  	<chromName><TAB><chromSize>

  	For example, Human (hg19):
  	chr1	249250621
  	chr2	243199373
  	...
  	chr18_gl000207_random	4262

  Tips: 
  	One can use the UCSC Genome Browser's MySQL database to extract
  	chromosome sizes. For example, H. sapiens:

  	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
  	"select chrom, size from hg19.chromInfo"  > hg19.genome


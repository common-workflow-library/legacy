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
  - complement
description: |
  Tool:    bedtools complement (aka complementBed)
  Version: v2.25.0
  Summary: Returns the base pair complement of a feature file.

  Usage:   bedtools complement [OPTIONS] -i <bed/gff/vcf> -g <genome>

  Notes: 
  	(1)  The genome file should tab delimited and structured as follows:
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


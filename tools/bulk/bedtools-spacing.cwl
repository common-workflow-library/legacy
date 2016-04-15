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
    description: '<bed/gff/vcf/bam>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#bed'
    type:
      - 'null'
      - boolean
    description: |
      If using BAM input, write output as BED.
    inputBinding:
      position: 1
      prefix: '-bed'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: |
      Print the header from the A file prior to results.
    inputBinding:
      position: 1
      prefix: '-header'
  - id: '#nobuf'
    type:
      - 'null'
      - boolean
    description: |
      Disable buffered output. Using this option will cause each line
      of output to be printed as it is generated, rather than saved
      in a buffer. This will make printing large output files
      noticeably slower, but can be useful in conjunction with
      other software tools and scripts that need to process one
      line of bedtools output at a time.
    inputBinding:
      position: 1
      prefix: '-nobuf'
  - id: '#iobuf'
    type:
      - 'null'
      - int
    description: |
      Specify amount of memory to use for input buffer.
      Takes an integer argument. Optional suffixes K/M/G supported.
      Note: currently has no effect with compressed files.
    inputBinding:
      position: 1
      prefix: '-iobuf'
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
  - spacing
description: |
  Tool:    bedtools spacing
  Version: v2.25.0
  Summary: Report (last col.) the gap lengths between intervals in a file.

  Usage:   bedtools spacing [OPTIONS] -i <bed/gff/vcf/bam>

  Notes: 
  	(1)  Input must be sorted by chrom,start (sort -k1,1 -k2,2n for BED).
  	(2)  The 1st element for each chrom will have NULL distance. (".").
  	(3)  Distance for overlapping intervaks is -1 and bookended is 0.

  Example: 
  	$ cat test.bed 
  	chr1    0   10 
  	chr1    10  20 
  	chr1    21  30 
  	chr1    35  45 
  	chr1    100 200 

  	$ bedtools spacing -i test.bed 
  	chr1    0   10  . 
  	chr1    10  20  0 
  	chr1    21  30  1 
  	chr1    35  45  5 
  	chr1    100 200 55 

  	-bed	If using BAM input, write output as BED.

  	-header	Print the header from the A file prior to results.

  	-nobuf	Disable buffered output. Using this option will cause each line
  		of output to be printed as it is generated, rather than saved
  		in a buffer. This will make printing large output files 
  		noticeably slower, but can be useful in conjunction with
  		other software tools and scripts that need to process one
  		line of bedtools output at a time.

  	-iobuf	Specify amount of memory to use for input buffer.
  		Takes an integer argument. Optional suffixes K/M/G supported.
  		Note: currently has no effect with compressed files.


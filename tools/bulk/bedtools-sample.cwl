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
  - id: '#n'
    type:
      - 'null'
      - int
    description: |
      The number of records to generate.
      - Default = 1,000,000.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#seed'
    type:
      - 'null'
      - int
    description: |
      Supply an integer seed for the shuffling.
      - By default, the seed is chosen automatically.
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: '-seed'
  - id: '#ubam'
    type:
      - 'null'
      - boolean
    description: |
      Write uncompressed BAM output. Default writes compressed BAM.
    inputBinding:
      position: 1
      prefix: '-ubam'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Require same strandedness.  That is, only give records
      that have the same strand. Use '-s forward' or '-s reverse'
      for forward or reverse strand records, respectively.
      - By default, records are reported without respect to strand.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: |
      Print the header from the input file prior to results.
    inputBinding:
      position: 1
      prefix: '-header'
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
  - sample
description: |
  Tool:    bedtools sample (aka sampleFile)
  Version: v2.25.0
  Summary: Take sample of input file(s) using reservoir sampling algorithm.

  Usage:   bedtools sample [OPTIONS] -i <bed/gff/vcf/bam>

  WARNING:	The current sample algorithm will hold all requested sample records in memory prior to output.
  		The user must ensure that there is adequate memory for this.

  Options: 
  	-n	The number of records to generate.
  		- Default = 1,000,000.
  		- (INTEGER)

  	-seed	Supply an integer seed for the shuffling.
  		- By default, the seed is chosen automatically.
  		- (INTEGER)

  	-ubam	Write uncompressed BAM output. Default writes compressed BAM.

  	-s	Require same strandedness.  That is, only give records
  		that have the same strand. Use '-s forward' or '-s reverse'
  		for forward or reverse strand records, respectively.
  		- By default, records are reported without respect to strand.

  	-header	Print the header from the input file prior to results.

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

  Notes:


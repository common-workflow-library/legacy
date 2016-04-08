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
  - id: '#hist'
    type:
      - 'null'
      - boolean
    description: |
      Report a histogram of coverage for each feature in A
      as well as a summary histogram for _all_ features in A.
      Output (tab delimited) after each feature in A:
      1) depth
      2) # bases at depth
      3) size of A
      4) % of A at depth
    inputBinding:
      position: 1
      prefix: '-hist'
  - id: '#d'
    type:
      - 'null'
      - boolean
    description: |
      Report the depth at each position in each A feature.
      Positions reported are one based.  Each position
      and depth follow the complete A feature.
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#counts'
    type:
      - 'null'
      - boolean
    description: |
      Only report the count of overlaps, don't compute fraction, etc.
    inputBinding:
      position: 1
      prefix: '-counts'
  - id: '#mean'
    type:
      - 'null'
      - boolean
    description: |
      Report the mean depth of all positions in each A feature.
    inputBinding:
      position: 1
      prefix: '-mean'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Require same strandedness.  That is, only report hits in B
      that overlap A on the _same_ strand.
      - By default, overlaps are reported without respect to strand.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: |
      Require different strandedness.  That is, only report hits in B
      that overlap A on the _opposite_ strand.
      - By default, overlaps are reported without respect to strand.
    inputBinding:
      position: 1
      prefix: '-S'
  - id: '#f'
    type:
      - 'null'
      - float
    description: |
      Minimum overlap required as a fraction of A.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#F'
    type:
      - 'null'
      - float
    description: |
      Minimum overlap required as a fraction of B.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)
    inputBinding:
      position: 1
      prefix: '-F'
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: |
      Require that the fraction overlap be reciprocal for A AND B.
      - In other words, if -f is 0.90 and -r is used, this requires
      that B overlap 90% of A and A _also_ overlaps 90% of B.
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#e'
    type:
      - 'null'
      - boolean
    description: |
      Require that the minimum fraction be satisfied for A OR B.
      - In other words, if -e is used with -f 0.90 and -F 0.10 this requires
      that either 90% of A is covered OR 10% of  B is covered.
      Without -e, both fractions would have to be satisfied.
    inputBinding:
      position: 1
      prefix: '-e'
  - id: '#split'
    type:
      - 'null'
      - boolean
    description: |
      Treat "split" BAM or BED12 entries as distinct BED intervals.
    inputBinding:
      position: 1
      prefix: '-split'
  - id: '#g'
    type:
      - 'null'
      - boolean
    description: |
      Provide a genome file to enforce consistent chromosome sort order
      across input files. Only applies when used with -sorted option.
    inputBinding:
      position: 1
      prefix: '-g'
  - id: '#nonamecheck'
    type:
      - 'null'
      - boolean
    description: >
      For sorted data, don't throw an error if the file has different naming
      conventions

      for the same chromosome. ex. "chr1" vs "chr01".
    inputBinding:
      position: 1
      prefix: '-nonamecheck'
  - id: '#sorted'
    type:
      - 'null'
      - boolean
    description: |
      Use the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input.
    inputBinding:
      position: 1
      prefix: '-sorted'
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
      Default Output:
      After each entry in A, reports:
      1) The number of features in B that overlapped the A interval.
      2) The number of bases in A that had non-zero coverage.
      3) The length of the entry in A.
      4) The fraction of bases in A that had non-zero coverage.
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
  - coverage
description: |
  Tool:    bedtools coverage (aka coverageBed)
  Version: v2.25.0
  Summary: Returns the depth and breadth of coverage of features from B
  	 on the intervals in A.

  Usage:   bedtools coverage [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

  Options: 
  	-hist	Report a histogram of coverage for each feature in A
  		as well as a summary histogram for _all_ features in A.

  		Output (tab delimited) after each feature in A:
  		  1) depth
  		  2) # bases at depth
  		  3) size of A
  		  4) % of A at depth

  	-d	Report the depth at each position in each A feature.
  		Positions reported are one based.  Each position
  		and depth follow the complete A feature.

  	-counts	Only report the count of overlaps, don't compute fraction, etc.

  	-mean	Report the mean depth of all positions in each A feature.

  	-s	Require same strandedness.  That is, only report hits in B
  		that overlap A on the _same_ strand.
  		- By default, overlaps are reported without respect to strand.

  	-S	Require different strandedness.  That is, only report hits in B
  		that overlap A on the _opposite_ strand.
  		- By default, overlaps are reported without respect to strand.

  	-f	Minimum overlap required as a fraction of A.
  		- Default is 1E-9 (i.e., 1bp).
  		- FLOAT (e.g. 0.50)

  	-F	Minimum overlap required as a fraction of B.
  		- Default is 1E-9 (i.e., 1bp).
  		- FLOAT (e.g. 0.50)

  	-r	Require that the fraction overlap be reciprocal for A AND B.
  		- In other words, if -f is 0.90 and -r is used, this requires
  		  that B overlap 90% of A and A _also_ overlaps 90% of B.

  	-e	Require that the minimum fraction be satisfied for A OR B.
  		- In other words, if -e is used with -f 0.90 and -F 0.10 this requires
  		  that either 90% of A is covered OR 10% of  B is covered.
  		  Without -e, both fractions would have to be satisfied.

  	-split	Treat "split" BAM or BED12 entries as distinct BED intervals.

  	-g	Provide a genome file to enforce consistent chromosome sort order
  		across input files. Only applies when used with -sorted option.

  	-nonamecheck	For sorted data, don't throw an error if the file has different naming conventions
  			for the same chromosome. ex. "chr1" vs "chr01".

  	-sorted	Use the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input.

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

  Default Output:  
  	 After each entry in A, reports: 
  	   1) The number of features in B that overlapped the A interval.
  	   2) The number of bases in A that had non-zero coverage.
  	   3) The length of the entry in A.
  	   4) The fraction of bases in A that had non-zero coverage.


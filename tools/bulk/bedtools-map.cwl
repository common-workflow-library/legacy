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
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: |
      Specify columns from the B file to map onto intervals in A.
      Default: 5.
      Multiple columns can be specified in a comma-delimited list.
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#o'
    type:
      - 'null'
      - boolean
    description: |
      Specify the operation that should be applied to -c.
      Valid operations:
      sum, min, max, absmin, absmax,
      mean, median,
      collapse (i.e., print a delimited list (duplicates allowed)),
      distinct (i.e., print a delimited list (NO duplicates allowed)),
      distinct_sort_num (as distinct, sorted numerically, ascending),
      distinct_sort_num_desc (as distinct, sorted numerically, desscending),
      distinct_only (delimited list of only unique values),
      count
      count_distinct (i.e., a count of the unique values in the column),
      first (i.e., just the first value in the column),
      last (i.e., just the last value in the column),
      Default: sum
      Multiple operations can be specified in a comma-delimited list.
      If there is only column, but multiple operations, all operations will be
      applied on that column. Likewise, if there is only one operation, but
      multiple columns, that operation will be applied to all columns.
      Otherwise, the number of columns must match the the number of operations,
      and will be applied in respective order.
      E.g., "-c 5,4,6 -o sum,mean,count" will give the sum of column 5,
      the mean of column 4, and the count of column 6.
      The order of output columns will match the ordering given in the command.
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#delim'
    type:
      - 'null'
      - boolean
    description: |
      Specify a custom delimiter for the collapse operations.
      - Example: -delim "|"
      - Default: ",".
    inputBinding:
      position: 1
      prefix: '-delim'
  - id: '#prec'
    type:
      - 'null'
      - boolean
    description: |
      Sets the decimal precision for output (Default: 5)
    inputBinding:
      position: 1
      prefix: '-prec'
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
  - map
description: |
  Tool:    bedtools map (aka mapBed)
  Version: v2.25.0
  Summary: Apply a function to a column from B intervals that overlap A.

  Usage:   bedtools map [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

  Options: 
  	-c	Specify columns from the B file to map onto intervals in A.
  		Default: 5.
  		Multiple columns can be specified in a comma-delimited list.

  	-o	Specify the operation that should be applied to -c.
  		Valid operations:
  		    sum, min, max, absmin, absmax,
  		    mean, median,
  		    collapse (i.e., print a delimited list (duplicates allowed)), 
  		    distinct (i.e., print a delimited list (NO duplicates allowed)), 
  		    distinct_sort_num (as distinct, sorted numerically, ascending),
  		    distinct_sort_num_desc (as distinct, sorted numerically, desscending),
  		    distinct_only (delimited list of only unique values),
  		    count
  		    count_distinct (i.e., a count of the unique values in the column), 
  		    first (i.e., just the first value in the column), 
  		    last (i.e., just the last value in the column), 
  		Default: sum
  		Multiple operations can be specified in a comma-delimited list.

  		If there is only column, but multiple operations, all operations will be
  		applied on that column. Likewise, if there is only one operation, but
  		multiple columns, that operation will be applied to all columns.
  		Otherwise, the number of columns must match the the number of operations,
  		and will be applied in respective order.
  		E.g., "-c 5,4,6 -o sum,mean,count" will give the sum of column 5,
  		the mean of column 4, and the count of column 6.
  		The order of output columns will match the ordering given in the command.


  	-delim	Specify a custom delimiter for the collapse operations.
  		- Example: -delim "|"
  		- Default: ",".

  	-prec	Sets the decimal precision for output (Default: 5)

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
  	(1) Both input files must be sorted by chrom, then start.


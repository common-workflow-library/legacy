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
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Force strandedness.  That is, only merge features
      that are on the same strand.
      - By default, merging is done without respect to strand.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: |
      Force merge for one specific strand only.
      Follow with + or - to force merge from only
      the forward or reverse strand, respectively.
      - By default, merging is done without respect to strand.
    inputBinding:
      position: 1
      prefix: '-S'
  - id: '#d'
    type:
      - 'null'
      - int
    description: |
      Maximum distance between features allowed for features
      to be merged.
      - Def. 0. That is, overlapping & book-ended features are merged.
      - (INTEGER)
      - Note: negative values enforce the number of b.p. required for overlap.
    inputBinding:
      position: 1
      prefix: '-d'
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
  - merge
description: |
  Tool:    bedtools merge (aka mergeBed)
  Version: v2.25.0
  Summary: Merges overlapping BED/GFF/VCF entries into a single interval.

  Usage:   bedtools merge [OPTIONS] -i <bed/gff/vcf>

  Options: 
  	-s	Force strandedness.  That is, only merge features
  		that are on the same strand.
  		- By default, merging is done without respect to strand.

  	-S	Force merge for one specific strand only.
  		Follow with + or - to force merge from only
  		the forward or reverse strand, respectively.
  		- By default, merging is done without respect to strand.

  	-d	Maximum distance between features allowed for features
  		to be merged.
  		- Def. 0. That is, overlapping & book-ended features are merged.
  		- (INTEGER)
  		- Note: negative values enforce the number of b.p. required for overlap.

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
  	(1) The input file (-i) file must be sorted by chrom, then start.


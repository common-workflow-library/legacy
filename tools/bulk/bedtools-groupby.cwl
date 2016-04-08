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
    type:
      - 'null'
      - boolean
    description: "\tInput file. Assumes \"stdin\" if omitted.\n"
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#g'
    type:
      - 'null'
      - boolean
    description: |
      -grp		Specify the columns (1-based) for the grouping.
      The columns must be comma separated.
      - Default: 1,2,3
    inputBinding:
      position: 1
      prefix: '-g'
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: "-opCols\tSpecify the column (1-based) that should be summarized.\n- Required.\n"
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#o'
    type:
      - 'null'
      - string
    description: |
      -ops		Specify the operation that should be applied to opCol.
      Valid operations:
      sum, count, count_distinct, min, max,
      mean, median, mode, antimode,
      stdev, sstdev (sample standard dev.),
      collapse (i.e., print a comma separated list (duplicates allowed)),
      distinct (i.e., print a comma separated list (NO duplicates allowed)),
      distinct_sort_num (as distinct, but sorted numerically, ascending),
      distinct_sort_num_desc (as distinct, but sorted numerically, descending),
      concat   (i.e., merge values into a single, non-delimited string),
      freqdesc (i.e., print desc. list of values:freq)
      freqasc (i.e., print asc. list of values:freq)
      first (i.e., print first value)
      last (i.e., print last value)
      - Default: sum
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
  - id: '#full'
    type:
      - 'null'
      - boolean
    description: "\tPrint all columns from input file.  The first line in the group is used.\nDefault: print only grouped columns.\n"
    inputBinding:
      position: 1
      prefix: '-full'
  - id: '#inheader'
    type:
      - 'null'
      - boolean
    description: |
      Input file has a header line - the first line will be ignored.
    inputBinding:
      position: 1
      prefix: '-inheader'
  - id: '#outheader'
    type:
      - 'null'
      - boolean
    description: |
      Print header line in the output, detailing the column names. 
      If the input file has headers (-inheader), the output file
      will use the input's column names.
      If the input file has no headers, the output file
      will use "col_1", "col_2", etc. as the column names.
    inputBinding:
      position: 1
      prefix: '-outheader'
  - id: '#header'
    type:
      - 'null'
      - boolean
    description: "\tsame as '-inheader -outheader'\n"
    inputBinding:
      position: 1
      prefix: '-header'
  - id: '#ignorecase'
    type:
      - 'null'
      - boolean
    description: |
      Group values regardless of upper/lower case.
    inputBinding:
      position: 1
      prefix: '-ignorecase'
  - id: '#prec'
    type:
      - 'null'
      - boolean
    description: |
      Sets the decimal precision for output (Default: 5)
    inputBinding:
      position: 1
      prefix: '-prec'
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
  - groupby
description: |
  Tool:    bedtools groupby 
  Version: v2.25.0
  Summary: Summarizes a dataset column based upon
  	 common column groupings. Akin to the SQL "group by" command.

  Usage:	 bedtools groupby -g [group_column(s)] -c [op_column(s)] -o [ops] 
  	 cat [FILE] | bedtools groupby -g [group_column(s)] -c [op_column(s)] -o [ops] 

  Options: 
  	-i		Input file. Assumes "stdin" if omitted.

  	-g -grp		Specify the columns (1-based) for the grouping.
  			The columns must be comma separated.
  			- Default: 1,2,3

  	-c -opCols	Specify the column (1-based) that should be summarized.
  			- Required.

  	-o -ops		Specify the operation that should be applied to opCol.
  			Valid operations:
  			    sum, count, count_distinct, min, max,
  			    mean, median, mode, antimode,
  			    stdev, sstdev (sample standard dev.),
  			    collapse (i.e., print a comma separated list (duplicates allowed)), 
  			    distinct (i.e., print a comma separated list (NO duplicates allowed)), 
  			    distinct_sort_num (as distinct, but sorted numerically, ascending), 
  			    distinct_sort_num_desc (as distinct, but sorted numerically, descending), 
  			    concat   (i.e., merge values into a single, non-delimited string), 
  			    freqdesc (i.e., print desc. list of values:freq)
  			    freqasc (i.e., print asc. list of values:freq)
  			    first (i.e., print first value)
  			    last (i.e., print last value)
  			- Default: sum

  		If there is only column, but multiple operations, all operations will be
  		applied on that column. Likewise, if there is only one operation, but
  		multiple columns, that operation will be applied to all columns.
  		Otherwise, the number of columns must match the the number of operations,
  		and will be applied in respective order.
  		E.g., "-c 5,4,6 -o sum,mean,count" will give the sum of column 5,
  		the mean of column 4, and the count of column 6.
  		The order of output columns will match the ordering given in the command.


  	-full		Print all columns from input file.  The first line in the group is used.
  			Default: print only grouped columns.

  	-inheader	Input file has a header line - the first line will be ignored.

  	-outheader	Print header line in the output, detailing the column names. 
  			If the input file has headers (-inheader), the output file
  			will use the input's column names.
  			If the input file has no headers, the output file
  			will use "col_1", "col_2", etc. as the column names.

  	-header		same as '-inheader -outheader'

  	-ignorecase	Group values regardless of upper/lower case.

  	-prec	Sets the decimal precision for output (Default: 5)

  	-delim	Specify a custom delimiter for the collapse operations.
  		- Example: -delim "|"
  		- Default: ",".

  Examples: 
  	$ cat ex1.out
  	chr1 10  20  A   chr1    15  25  B.1 1000    ATAT
  	chr1 10  20  A   chr1    25  35  B.2 10000   CGCG

  	$ groupBy -i ex1.out -g 1,2,3,4 -c 9 -o sum
  	chr1 10  20  A   11000

  	$ groupBy -i ex1.out -grp 1,2,3,4 -opCols 9,9 -ops sum,max
  	chr1 10  20  A   11000   10000

  	$ groupBy -i ex1.out -g 1,2,3,4 -c 8,9 -o collapse,mean
  	chr1 10  20  A   B.1,B.2,    5500

  	$ cat ex1.out | groupBy -g 1,2,3,4 -c 8,9 -o collapse,mean
  	chr1 10  20  A   B.1,B.2,    5500

  	$ cat ex1.out | groupBy -g 1,2,3,4 -c 10 -o concat
  	chr1 10  20  A   ATATCGCG

  Notes: 
  	(1)  The input file/stream should be sorted/grouped by the -grp. columns
  	(2)  If -i is unspecified, input is assumed to come from stdin.


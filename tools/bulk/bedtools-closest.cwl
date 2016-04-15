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
  - id: '#d'
    type:
      - 'null'
      - boolean
    description: |
      In addition to the closest feature in B, 
      report its distance to A as an extra column.
      - The reported distance for overlapping features will be 0.
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#D'
    type:
      - 'null'
      - boolean
    description: |
      Like -d, report the closest feature in B, and its distance to A
      as an extra column. Unlike -d, use negative distances to report
      upstream features.
      The options for defining which orientation is "upstream" are:
      - "ref"   Report distance with respect to the reference genome.
      B features with a lower (start, stop) are upstream
      - "a"     Report distance with respect to A.
      When A is on the - strand, "upstream" means B has a
      higher (start,stop).
      - "b"     Report distance with respect to B.
      When B is on the - strand, "upstream" means A has a
      higher (start,stop).
    inputBinding:
      position: 1
      prefix: '-D'
  - id: '#io'
    type:
      - 'null'
      - boolean
    description: |
      Ignore features in B that overlap A.  That is, we want close,
      yet not touching features only.
    inputBinding:
      position: 1
      prefix: '-io'
  - id: '#iu'
    type:
      - 'null'
      - boolean
    description: |
      Ignore features in B that are upstream of features in A.
      This option requires -D and follows its orientation
      rules for determining what is "upstream".
    inputBinding:
      position: 1
      prefix: '-iu'
  - id: '#id'
    type:
      - 'null'
      - boolean
    description: |
      Ignore features in B that are downstream of features in A.
      This option requires -D and follows its orientation
      rules for determining what is "downstream".
    inputBinding:
      position: 1
      prefix: '-id'
  - id: '#fu'
    type:
      - 'null'
      - boolean
    description: |
      Choose first from features in B that are upstream of features in A.
      This option requires -D and follows its orientation
      rules for determining what is "upstream".
    inputBinding:
      position: 1
      prefix: '-fu'
  - id: '#fd'
    type:
      - 'null'
      - boolean
    description: |
      Choose first from features in B that are downstream of features in A.
      This option requires -D and follows its orientation
      rules for determining what is "downstream".
    inputBinding:
      position: 1
      prefix: '-fd'
  - id: '#t'
    type:
      - 'null'
      - boolean
    description: |
      How ties for closest feature are handled.  This occurs when two
      features in B have exactly the same "closeness" with A.
      By default, all such features in B are reported.
      Here are all the options:
      - "all"    Report all ties (default).
      - "first"  Report the first tie that occurred in the B file.
      - "last"   Report the last tie that occurred in the B file.
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#mdb'
    type:
      - 'null'
      - boolean
    description: |
      How multiple databases are resolved.
      - "each"    Report closest records for each database (default).
      - "all"  Report closest records among all databases.
    inputBinding:
      position: 1
      prefix: '-mdb'
  - id: '#k'
    type:
      - 'null'
      - boolean
    description: "Report the k closest hits. Default is 1. If tieMode = \"all\", \n- all ties will still be reported.\n"
    inputBinding:
      position: 1
      prefix: '-k'
  - id: '#N'
    type:
      - 'null'
      - boolean
    description: |
      Require that the query and the closest hit have different names.
      For BED, the 4th column is compared.
    inputBinding:
      position: 1
      prefix: '-N'
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
  - id: '#names'
    type:
      - 'null'
      - boolean
    description: |
      When using multiple databases, provide an alias for each that
      will appear instead of a fileId when also printing the DB record.
    inputBinding:
      position: 1
      prefix: '-names'
  - id: '#filenames'
    type:
      - 'null'
      - boolean
    description: |
      When using multiple databases, show each complete filename
      instead of a fileId when also printing the DB record.
    inputBinding:
      position: 1
      prefix: '-filenames'
  - id: '#sortout'
    type:
      - 'null'
      - boolean
    description: |
      When using multiple databases, sort the output DB hits
      for each record.
    inputBinding:
      position: 1
      prefix: '-sortout'
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
  - closest
description: |
  Tool:    bedtools closest (aka closestBed)
  Version: v2.25.0
  Summary: For each feature in A, finds the closest 
  	 feature (upstream or downstream) in B.

  Usage:   bedtools closest [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

  Options: 
  	-d	In addition to the closest feature in B, 
  		report its distance to A as an extra column.
  		- The reported distance for overlapping features will be 0.

  	-D	Like -d, report the closest feature in B, and its distance to A
  		as an extra column. Unlike -d, use negative distances to report
  		upstream features.
  		The options for defining which orientation is "upstream" are:
  		- "ref"   Report distance with respect to the reference genome. 
  		            B features with a lower (start, stop) are upstream
  		- "a"     Report distance with respect to A.
  		            When A is on the - strand, "upstream" means B has a
  		            higher (start,stop).
  		- "b"     Report distance with respect to B.
  		            When B is on the - strand, "upstream" means A has a
  		            higher (start,stop).

  	-io	Ignore features in B that overlap A.  That is, we want close,
  		yet not touching features only.

  	-iu	Ignore features in B that are upstream of features in A.
  		This option requires -D and follows its orientation
  		rules for determining what is "upstream".

  	-id	Ignore features in B that are downstream of features in A.
  		This option requires -D and follows its orientation
  		rules for determining what is "downstream".

  	-fu	Choose first from features in B that are upstream of features in A.
  		This option requires -D and follows its orientation
  		rules for determining what is "upstream".

  	-fd	Choose first from features in B that are downstream of features in A.
  		This option requires -D and follows its orientation
  		rules for determining what is "downstream".

  	-t	How ties for closest feature are handled.  This occurs when two
  		features in B have exactly the same "closeness" with A.
  		By default, all such features in B are reported.
  		Here are all the options:
  		- "all"    Report all ties (default).
  		- "first"  Report the first tie that occurred in the B file.
  		- "last"   Report the last tie that occurred in the B file.

  	-mdb	How multiple databases are resolved.
  		- "each"    Report closest records for each database (default).
  		- "all"  Report closest records among all databases.

  	-k	Report the k closest hits. Default is 1. If tieMode = "all", 
  		- all ties will still be reported.

  	-N	Require that the query and the closest hit have different names.
  		For BED, the 4th column is compared.

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

  	-names	When using multiple databases, provide an alias for each that
  		will appear instead of a fileId when also printing the DB record.

  	-filenames	When using multiple databases, show each complete filename
  			instead of a fileId when also printing the DB record.

  	-sortout	When using multiple databases, sort the output DB hits
  			for each record.

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
  	Reports "none" for chrom and "-1" for all other fields when a feature
  	is not found in B on the same chromosome as the feature in A.
  	E.g. none	-1	-1


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
    description: '<bed/gff/vcf/bam>'
    inputBinding:
      position: 3
      prefix: '-b'
  - id: '#a'
    type: File
    description: '<bed/gff/vcf/bam>'
    inputBinding:
      position: 2
      prefix: '-a'
  - id: '#wa'
    type:
      - 'null'
      - boolean
    description: |
      Write the original entry in A for each overlap.
    inputBinding:
      position: 1
      prefix: '-wa'
  - id: '#wb'
    type:
      - 'null'
      - boolean
    description: |
      Write the original entry in B for each overlap.
      - Useful for knowing _what_ A overlaps. Restricted by -f and -r.
    inputBinding:
      position: 1
      prefix: '-wb'
  - id: '#loj'
    type:
      - 'null'
      - boolean
    description: |
      Perform a "left outer join". That is, for each feature in A
      report each overlap with B.  If no overlaps are found,
      report a NULL feature for B.
    inputBinding:
      position: 1
      prefix: '-loj'
  - id: '#wo'
    type:
      - 'null'
      - boolean
    description: |
      Write the original A and B entries plus the number of base
      pairs of overlap between the two features.
      - Overlaps restricted by -f and -r.
      Only A features with overlap are reported.
    inputBinding:
      position: 1
      prefix: '-wo'
  - id: '#wao'
    type:
      - 'null'
      - boolean
    description: |
      Write the original A and B entries plus the number of base
      pairs of overlap between the two features.
      - Overlapping features restricted by -f and -r.
      However, A features w/o overlap are also reported
      with a NULL B feature and overlap = 0.
    inputBinding:
      position: 1
      prefix: '-wao'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: |
      Write the original A entry _once_ if _any_ overlaps found in B.
      - In other words, just report the fact >=1 hit was found.
      - Overlaps restricted by -f and -r.
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
      - Overlaps restricted by -f and -r.
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#v'
    type:
      - 'null'
      - boolean
    description: |
      Only report those entries in A that have _no overlaps_ with B.
      - Similar to "grep -v" (an homage).
    inputBinding:
      position: 1
      prefix: '-v'
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
  - intersect
description: |
  Tool:    bedtools intersect (aka intersectBed)
  Version: v2.25.0
  Summary: Report overlaps between two feature files.

  Usage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>

  	Note: -b may be followed with multiple databases and/or 
  	wildcard (*) character(s). 
  Options: 
  	-wa	Write the original entry in A for each overlap.

  	-wb	Write the original entry in B for each overlap.
  		- Useful for knowing _what_ A overlaps. Restricted by -f and -r.

  	-loj	Perform a "left outer join". That is, for each feature in A
  		report each overlap with B.  If no overlaps are found, 
  		report a NULL feature for B.

  	-wo	Write the original A and B entries plus the number of base
  		pairs of overlap between the two features.
  		- Overlaps restricted by -f and -r.
  		  Only A features with overlap are reported.

  	-wao	Write the original A and B entries plus the number of base
  		pairs of overlap between the two features.
  		- Overlapping features restricted by -f and -r.
  		  However, A features w/o overlap are also reported
  		  with a NULL B feature and overlap = 0.

  	-u	Write the original A entry _once_ if _any_ overlaps found in B.
  		- In other words, just report the fact >=1 hit was found.
  		- Overlaps restricted by -f and -r.

  	-c	For each entry in A, report the number of overlaps with B.
  		- Reports 0 for A entries that have no overlap with B.
  		- Overlaps restricted by -f and -r.

  	-v	Only report those entries in A that have _no overlaps_ with B.
  		- Similar to "grep -v" (an homage).

  	-ubam	Write uncompressed BAM output. Default writes compressed BAM.

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
  	(1) When a BAM file is used for the A file, the alignment is retained if overlaps exist,
  	and exlcuded if an overlap cannot be found.  If multiple overlaps exist, they are not
  	reported, as we are only testing for one or more overlaps.


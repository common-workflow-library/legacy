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
  - id: '#alnnbam'
    type: File
    description: aln.n.bam
    inputBinding:
      position: 5
  - id: '#'
    type: boolean
    description: ...
    inputBinding:
      position: 4
  - id: '#aln2bam'
    type: File
    description: aln.2.bam
    inputBinding:
      position: 3
  - id: '#bams'
    type:
      - 'null'
      - boolean
    description: |
      The bam files.
    inputBinding:
      position: 1
      prefix: '-bams'
  - id: '#bed'
    type:
      - 'null'
      - boolean
    description: |
      The bed file.
    inputBinding:
      position: 1
      prefix: '-bed'
  - id: '#split'
    type:
      - 'null'
      - boolean
    description: |
      Treat "split" BAM or BED12 entries as distinct BED intervals.
    inputBinding:
      position: 1
      prefix: '-split'
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
      Minimum overlap required as a fraction of each A.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: |
      Require that the fraction overlap be reciprocal for each A and B.
      - In other words, if -f is 0.90 and -r is used, this requires
      that B overlap 90% of A and A _also_ overlaps 90% of B.
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#q'
    type:
      - 'null'
      - boolean
    description: |
      Minimum mapping quality allowed. Default is 0.
    inputBinding:
      position: 1
      prefix: '-q'
  - id: '#D'
    type:
      - 'null'
      - boolean
    description: |
      Include duplicate reads.  Default counts non-duplicates only
    inputBinding:
      position: 1
      prefix: '-D'
  - id: '#F'
    type:
      - 'null'
      - boolean
    description: |
      Include failed-QC reads.  Default counts pass-QC reads only
    inputBinding:
      position: 1
      prefix: '-F'
  - id: '#p'
    type:
      - 'null'
      - boolean
    description: |
      Only count proper pairs.  Default counts all alignments with
      MAPQ > -q argument, regardless of the BAM FLAG field.
    inputBinding:
      position: 1
      prefix: '-p'
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
  - multicov
description: |
  Tool:    bedtools multicov (aka multiBamCov)
  Version: v2.25.0
  Summary: Counts sequence coverage for multiple bams at specific loci.

  Usage:   bedtools multicov [OPTIONS] -bams aln.1.bam aln.2.bam ... aln.n.bam -bed <bed/gff/vcf>

  Options: 
  	-bams	The bam files.

  	-bed	The bed file.

  	-split	Treat "split" BAM or BED12 entries as distinct BED intervals.

  	-s	Require same strandedness.  That is, only report hits in B
  		that overlap A on the _same_ strand.
  		- By default, overlaps are reported without respect to strand.

  	-S	Require different strandedness.  That is, only report hits in B
  		that overlap A on the _opposite_ strand.
  		- By default, overlaps are reported without respect to strand.

  	-f	Minimum overlap required as a fraction of each A.
  		- Default is 1E-9 (i.e., 1bp).
  		- FLOAT (e.g. 0.50)

  	-r	Require that the fraction overlap be reciprocal for each A and B.
  		- In other words, if -f is 0.90 and -r is used, this requires
  		  that B overlap 90% of A and A _also_ overlaps 90% of B.

  	-q	Minimum mapping quality allowed. Default is 0.

  	-D	Include duplicate reads.  Default counts non-duplicates only

  	-F	Include failed-QC reads.  Default counts pass-QC reads only

  	-p	Only count proper pairs.  Default counts all alignments with
  		MAPQ > -q argument, regardless of the BAM FLAG field.


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
    description: '<bedpe>'
    inputBinding:
      position: 2
      prefix: '-a'
  - id: '#abam'
    type:
      - 'null'
      - boolean
    description: >
      The A input file is in BAM format.  Output will be BAM as well. Replaces
      -a.

      - Requires BAM to be grouped or sorted by query.
    inputBinding:
      position: 1
      prefix: '-abam'
  - id: '#ubam'
    type:
      - 'null'
      - boolean
    description: |
      Write uncompressed BAM output. Default writes compressed BAM.
      is to write output in BAM when using -abam.
    inputBinding:
      position: 1
      prefix: '-ubam'
  - id: '#bedpe'
    type:
      - 'null'
      - boolean
    description: |
      When using BAM input (-abam), write output as BEDPE. The default
      is to write output in BAM when using -abam.
    inputBinding:
      position: 1
      prefix: '-bedpe'
  - id: '#ed'
    type:
      - 'null'
      - boolean
    description: |
      Use BAM total edit distance (NM tag) for BEDPE score.
      - Default for BEDPE is to use the minimum of
      of the two mapping qualities for the pair.
      - When -ed is used the total edit distance
      from the two mates is reported as the score.
    inputBinding:
      position: 1
      prefix: '-ed'
  - id: '#f'
    type:
      - 'null'
      - boolean
    description: |
      Minimum overlap required as fraction of A (e.g. 0.05).
      Default is 1E-9 (effectively 1bp).
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Require same strandedness when finding overlaps.
      Default is to ignore stand.
      Not applicable with -type inspan or -type outspan.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: |
      Require different strandedness when finding overlaps.
      Default is to ignore stand.
      Not applicable with -type inspan or -type outspan.
    inputBinding:
      position: 1
      prefix: '-S'
  - id: '#type'
    type:
      - 'null'
      - boolean
    description: |
      	Approach to reporting overlaps between BEDPE and BED.
      either	Report overlaps if either end of A overlaps B.
      - Default.
      neither	Report A if neither end of A overlaps B.
      both	Report overlaps if both ends of A overlap  B.
      xor	Report overlaps if one and only one end of A overlaps B.
      notboth	Report overlaps if neither end or one and only one
      end of A overlap B.  That is, xor + neither.
      ispan	Report overlaps between [end1, start2] of A and B.
      - Note: If chrom1 <> chrom2, entry is ignored.
      ospan	Report overlaps between [start1, end2] of A and B.
      - Note: If chrom1 <> chrom2, entry is ignored.
      notispan	Report A if ispan of A doesn't overlap B.
      - Note: If chrom1 <> chrom2, entry is ignored.
      notospan	Report A if ospan of A doesn't overlap B.
      - Note: If chrom1 <> chrom2, entry is ignored.
      Refer to the BEDTools manual for BEDPE format.
    inputBinding:
      position: 1
      prefix: '-type'
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
  - pairtobed
description: |
  Tool:    bedtools pairtobed (aka pairToBed)
  Version: v2.25.0
  Summary: Report overlaps between a BEDPE file and a BED/GFF/VCF file.

  Usage:   bedtools pairtobed [OPTIONS] -a <bedpe> -b <bed/gff/vcf>

  Options: 
  	-abam	The A input file is in BAM format.  Output will be BAM as well. Replaces -a.
  		- Requires BAM to be grouped or sorted by query.

  	-ubam	Write uncompressed BAM output. Default writes compressed BAM.

  		is to write output in BAM when using -abam.

  	-bedpe	When using BAM input (-abam), write output as BEDPE. The default
  		is to write output in BAM when using -abam.

  	-ed	Use BAM total edit distance (NM tag) for BEDPE score.
  		- Default for BEDPE is to use the minimum of
  		  of the two mapping qualities for the pair.
  		- When -ed is used the total edit distance
  		  from the two mates is reported as the score.

  	-f	Minimum overlap required as fraction of A (e.g. 0.05).
  		Default is 1E-9 (effectively 1bp).

  	-s	Require same strandedness when finding overlaps.
  		Default is to ignore stand.
  		Not applicable with -type inspan or -type outspan.

  	-S	Require different strandedness when finding overlaps.
  		Default is to ignore stand.
  		Not applicable with -type inspan or -type outspan.

  	-type 	Approach to reporting overlaps between BEDPE and BED.

  		either	Report overlaps if either end of A overlaps B.
  			- Default.
  		neither	Report A if neither end of A overlaps B.
  		both	Report overlaps if both ends of A overlap  B.
  		xor	Report overlaps if one and only one end of A overlaps B.
  		notboth	Report overlaps if neither end or one and only one 
  			end of A overlap B.  That is, xor + neither.

  		ispan	Report overlaps between [end1, start2] of A and B.
  			- Note: If chrom1 <> chrom2, entry is ignored.

  		ospan	Report overlaps between [start1, end2] of A and B.
  			- Note: If chrom1 <> chrom2, entry is ignored.

  		notispan	Report A if ispan of A doesn't overlap B.
  				- Note: If chrom1 <> chrom2, entry is ignored.

  		notospan	Report A if ospan of A doesn't overlap B.
  				- Note: If chrom1 <> chrom2, entry is ignored.

  Refer to the BEDTools manual for BEDPE format.


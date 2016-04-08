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
  - id: '#g'
    type: boolean
    description: '<genome>'
    inputBinding:
      position: 3
      prefix: '-g'
  - id: '#i'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#excl'
    type:
      - 'null'
      - boolean
    description: |
      A BED/GFF/VCF file of coordinates in which features in -i
      should not be placed (e.g. gaps.bed).
    inputBinding:
      position: 1
      prefix: '-excl'
  - id: '#incl'
    type:
      - 'null'
      - boolean
    description: |
      Instead of randomly placing features in a genome, the -incl
      options defines a BED/GFF/VCF file of coordinates in which
      features in -i should be randomly placed (e.g. genes.bed).
    inputBinding:
      position: 1
      prefix: '-incl'
  - id: '#chrom'
    type:
      - 'null'
      - boolean
    description: |
      Keep features in -i on the same chromosome.
      - By default, the chrom and position are randomly chosen.
      - NOTE: Forces use of -chromFirst (see below).
    inputBinding:
      position: 1
      prefix: '-chrom'
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
  - id: '#f'
    type:
      - 'null'
      - float
    description: |
      Maximum overlap (as a fraction of the -i feature) with an -excl
      feature that is tolerated before searching for a new,
      randomized locus. For example, -f 0.10 allows up to 10%
      of a randomized feature to overlap with a given feature
      in the -excl file. **Cannot be used with -incl file.**
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#chromFirst'
    type:
      - 'null'
      - boolean
    description: |
      Instead of choosing a position randomly among the entire
      genome (the default), first choose a chrom randomly, and then
      choose a random start coordinate on that chrom.  This leads
      to features being ~uniformly distributed among the chroms,
      as opposed to features being distribute as a function of chrom size.
    inputBinding:
      position: 1
      prefix: '-chromFirst'
  - id: '#bedpe'
    type:
      - 'null'
      - boolean
    description: |
      Indicate that the A file is in BEDPE format.
    inputBinding:
      position: 1
      prefix: '-bedpe'
  - id: '#maxTries'
    type:
      - 'null'
      - boolean
    description: |
      Max. number of attempts to find a home for a shuffled interval
      in the presence of -incl or -excl.
      Default = 1000.
    inputBinding:
      position: 1
      prefix: '-maxTries'
  - id: '#noOverlapping'
    type:
      - 'null'
      - boolean
    description: |
      Don't allow shuffled intervals to overlap.
    inputBinding:
      position: 1
      prefix: '-noOverlapping'
  - id: '#allowBeyondChromEnd'
    type:
      - 'null'
      - boolean
    description: |
      Allow shuffled intervals to be relocated to a position
      in which the entire original interval cannot fit w/o exceeding
      the end of the chromosome.  In this case, the end coordinate of the
      shuffled interval will be set to the chromosome's length.
      By default, an interval's original length must be fully-contained
      within the chromosome.
    inputBinding:
      position: 1
      prefix: '-allowBeyondChromEnd'
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
  - shuffle
description: |
  Tool:    bedtools shuffle (aka shuffleBed)
  Version: v2.25.0
  Summary: Randomly permute the locations of a feature file among a genome.

  Usage:   bedtools shuffle [OPTIONS] -i <bed/gff/vcf> -g <genome>

  Options: 
  	-excl	A BED/GFF/VCF file of coordinates in which features in -i
  		should not be placed (e.g. gaps.bed).

  	-incl	Instead of randomly placing features in a genome, the -incl
  		options defines a BED/GFF/VCF file of coordinates in which 
  		features in -i should be randomly placed (e.g. genes.bed). 
  	-chrom	Keep features in -i on the same chromosome.
  		- By default, the chrom and position are randomly chosen.
  		- NOTE: Forces use of -chromFirst (see below).

  	-seed	Supply an integer seed for the shuffling.
  		- By default, the seed is chosen automatically.
  		- (INTEGER)

  	-f	Maximum overlap (as a fraction of the -i feature) with an -excl
  		feature that is tolerated before searching for a new, 
  		randomized locus. For example, -f 0.10 allows up to 10%
  		of a randomized feature to overlap with a given feature
  		in the -excl file. **Cannot be used with -incl file.**
  		- Default is 1E-9 (i.e., 1bp).
  		- FLOAT (e.g. 0.50)

  	-chromFirst	
  		Instead of choosing a position randomly among the entire
  		genome (the default), first choose a chrom randomly, and then
  		choose a random start coordinate on that chrom.  This leads
  		to features being ~uniformly distributed among the chroms,
  		as opposed to features being distribute as a function of chrom size.

  	-bedpe	Indicate that the A file is in BEDPE format.

  	-maxTries	
  		Max. number of attempts to find a home for a shuffled interval
  		in the presence of -incl or -excl.
  		Default = 1000.
  	-noOverlapping	
  		Don't allow shuffled intervals to overlap.
  	-allowBeyondChromEnd	
  		Allow shuffled intervals to be relocated to a position
  		in which the entire original interval cannot fit w/o exceeding
  		the end of the chromosome.  In this case, the end coordinate of the
  		shuffled interval will be set to the chromosome's length.
  		By default, an interval's original length must be fully-contained
  		within the chromosome.
  Notes: 
  	(1)  The genome file should tab delimited and structured as follows:
  	     <chromName><TAB><chromSize>

  	For example, Human (hg19):
  	chr1	249250621
  	chr2	243199373
  	...
  	chr18_gl000207_random	4262

  Tips: 
  	One can use the UCSC Genome Browser's MySQL database to extract
  	chromosome sizes. For example, H. sapiens:

  	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
  	"select chrom, size from hg19.chromInfo"  > hg19.genome


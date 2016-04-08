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
  - id: '#fi'
    type:
      - 'null'
      - boolean
    description: |
      Input FASTA file
    inputBinding:
      position: 1
      prefix: '-fi'
  - id: '#bed'
    type:
      - 'null'
      - boolean
    description: |
      BED/GFF/VCF file of ranges to extract from -fi
    inputBinding:
      position: 1
      prefix: '-bed'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      Profile the sequence according to strand.
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#seq'
    type:
      - 'null'
      - boolean
    description: |
      Print the extracted sequence
    inputBinding:
      position: 1
      prefix: '-seq'
  - id: '#pattern'
    type:
      - 'null'
      - boolean
    description: |
      Report the number of times a user-defined sequence
      is observed (case-sensitive).
    inputBinding:
      position: 1
      prefix: '-pattern'
  - id: '#C'
    type:
      - 'null'
      - boolean
    description: |
      Ignore case when matching -pattern. By defaulty, case matters.
    inputBinding:
      position: 1
      prefix: '-C'
  - id: '#fullHeader'
    type:
      - 'null'
      - boolean
    description: |
      Use full fasta header.
      - By default, only the word before the first space or tab is used.
      Output format:
      The following information will be reported after each BED entry:
      1) %AT content
      2) %GC content
      3) Number of As observed
      4) Number of Cs observed
      5) Number of Gs observed
      6) Number of Ts observed
      7) Number of Ns observed
      8) Number of other bases observed
      9) The length of the explored sequence/interval.
      10) The seq. extracted from the FASTA file. (opt., if -seq is used)
      11) The number of times a user's pattern was observed.
      (opt., if -pattern is used.)
    inputBinding:
      position: 1
      prefix: '-fullHeader'
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
  - nuc
description: |
  Tool:    bedtools nuc (aka nucBed)
  Version: v2.25.0
  Summary: Profiles the nucleotide content of intervals in a fasta file.

  Usage:   bedtools nuc [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>

  Options: 
  	-fi	Input FASTA file

  	-bed	BED/GFF/VCF file of ranges to extract from -fi

  	-s	Profile the sequence according to strand.

  	-seq	Print the extracted sequence

  	-pattern	Report the number of times a user-defined sequence
  			is observed (case-sensitive).

  	-C	Ignore case when matching -pattern. By defaulty, case matters.

  	-fullHeader	Use full fasta header.
  		- By default, only the word before the first space or tab is used.

  Output format: 
  	The following information will be reported after each BED entry:
  	    1) %AT content
  	    2) %GC content
  	    3) Number of As observed
  	    4) Number of Cs observed
  	    5) Number of Gs observed
  	    6) Number of Ts observed
  	    7) Number of Ns observed
  	    8) Number of other bases observed
  	    9) The length of the explored sequence/interval.
  	    10) The seq. extracted from the FASTA file. (opt., if -seq is used)
  	    11) The number of times a user's pattern was observed.
  	        (opt., if -pattern is used.)


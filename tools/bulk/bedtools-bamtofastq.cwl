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
  - id: '#fq'
    type: boolean
    description: '<FQ>'
    inputBinding:
      position: 4
      prefix: '-fq'
  - id: '#BAM'
    type: File
    description: '<BAM>'
    inputBinding:
      position: 3
  - id: '#i'
    type: boolean
    description: ''
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#fq2'
    type:
      - 'null'
      - boolean
    description: |
      FASTQ for second end.  Used if BAM contains paired-end data.
      BAM should be sorted by query name is creating paired FASTQ.
    inputBinding:
      position: 1
      prefix: '-fq2'
  - id: '#tags'
    type:
      - 'null'
      - boolean
    description: |
      Create FASTQ based on the mate info
      in the BAM R2 and Q2 tags.
    inputBinding:
      position: 1
      prefix: '-tags'
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
  - bamtofastq
description: |
  Tool:    bedtools bamtofastq (aka bamToFastq)
  Version: v2.25.0
  Summary: Convert BAM alignments to FASTQ files. 

  Usage:   bamToFastq [OPTIONS] -i <BAM> -fq <FQ> 

  Options:
  	-fq2	FASTQ for second end.  Used if BAM contains paired-end data.
  		BAM should be sorted by query name is creating paired FASTQ.

  	-tags	Create FASTQ based on the mate info
  		in the BAM R2 and Q2 tags.

  Tips: 
  	If you want to create a single, interleaved FASTQ file 
  	for paired-end data, you can just write both to /dev/stdout:

  	bedtools bamtofastq -i x.bam -fq /dev/stdout -fq2 /dev/stdout > x.ilv.fq


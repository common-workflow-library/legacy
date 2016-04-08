#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: bwa-docker.cwl
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#in2fq'
    type: File
    description: '<in2.fq>'
    inputBinding:
      position: 6
  - id: '#in1fq'
    type: File
    description: '<in1.fq>'
    inputBinding:
      position: 5
  - id: '#in2sai'
    type: File
    description: '<in2.sai>'
    inputBinding:
      position: 4
  - id: '#in1sai'
    type: File
    description: '<in1.sai>'
    inputBinding:
      position: 3
  - id: '#prefix'
    type: boolean
    description: '<prefix>'
    inputBinding:
      position: 2
  - id: '#a'
    type:
      - 'null'
      - int
    description: |
      INT   maximum insert size [500]
    inputBinding:
      position: 1
      prefix: '-a'
  - id: '#o'
    type:
      - 'null'
      - int
    description: |
      INT   maximum occurrences for one end [100000]
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#n'
    type:
      - 'null'
      - int
    description: |
      INT   maximum hits to output for paired reads [3]
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#N'
    type:
      - 'null'
      - int
    description: |
      INT   maximum hits to output for discordant pairs [10]
    inputBinding:
      position: 1
      prefix: '-N'
  - id: '#c'
    type:
      - 'null'
      - float
    description: |
      FLOAT prior of chimeric rate (lower bound) [1.0e-05]
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#f'
    type:
      - 'null'
      - File
    description: |
      FILE  sam file to output results to [stdout]
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#r'
    type:
      - 'null'
      - string
    description: "STR   read group header line such as `@RG\\tID:foo\\tSM:bar' [null]\n"
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#P'
    type:
      - 'null'
      - boolean
    description: "      preload index into memory (for base-space reads only)\n"
    inputBinding:
      position: 1
      prefix: '-P'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: "      disable Smith-Waterman for the unmapped mate\n"
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#A'
    type:
      - 'null'
      - boolean
    description: |
            disable insert size estimate (force -s)
      Notes: 1. For SOLiD reads, <in1.fq> corresponds R3 reads and <in2.fq> to F3.
      2. For reads shorter than 30bp, applying a smaller -o is recommended to
      to get a sensible speed at the cost of pairing accuracy.
    inputBinding:
      position: 1
      prefix: '-A'
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
  - bwa
  - sampe
description: |
  Usage:   bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq>

  Options: -a INT   maximum insert size [500]
           -o INT   maximum occurrences for one end [100000]
           -n INT   maximum hits to output for paired reads [3]
           -N INT   maximum hits to output for discordant pairs [10]
           -c FLOAT prior of chimeric rate (lower bound) [1.0e-05]
           -f FILE  sam file to output results to [stdout]
           -r STR   read group header line such as `@RG\tID:foo\tSM:bar' [null]
           -P       preload index into memory (for base-space reads only)
           -s       disable Smith-Waterman for the unmapped mate
           -A       disable insert size estimate (force -s)

  Notes: 1. For SOLiD reads, <in1.fq> corresponds R3 reads and <in2.fq> to F3.
         2. For reads shorter than 30bp, applying a smaller -o is recommended to
            to get a sensible speed at the cost of pairing accuracy.


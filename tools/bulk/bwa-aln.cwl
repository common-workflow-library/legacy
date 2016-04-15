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
  - id: '#infq'
    type: File
    description: '<in.fq>'
    inputBinding:
      position: 3
  - id: '#prefix'
    type: boolean
    description: '<prefix>'
    inputBinding:
      position: 2
  - id: '#n'
    type:
      - 'null'
      - float
    description: |
      NUM    max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#o'
    type:
      - 'null'
      - int
    description: |
      INT    maximum number or fraction of gap opens [1]
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#e'
    type:
      - 'null'
      - int
    description: |
      INT    maximum number of gap extensions, -1 for disabling long gaps [-1]
    inputBinding:
      position: 1
      prefix: '-e'
  - id: '#i'
    type:
      - 'null'
      - int
    description: |
      INT    do not put an indel within INT bp towards the ends [5]
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#d'
    type:
      - 'null'
      - int
    description: |
      INT    maximum occurrences for extending a long deletion [10]
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      INT    seed length [32]
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#k'
    type:
      - 'null'
      - int
    description: |
      INT    maximum differences in the seed [2]
    inputBinding:
      position: 1
      prefix: '-k'
  - id: '#m'
    type:
      - 'null'
      - int
    description: |
      INT    maximum entries in the queue [2000000]
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#t'
    type:
      - 'null'
      - int
    description: |
      INT    number of threads [1]
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#M'
    type:
      - 'null'
      - int
    description: |
      INT    mismatch penalty [3]
    inputBinding:
      position: 1
      prefix: '-M'
  - id: '#O'
    type:
      - 'null'
      - int
    description: |
      INT    gap open penalty [11]
    inputBinding:
      position: 1
      prefix: '-O'
  - id: '#E'
    type:
      - 'null'
      - int
    description: |
      INT    gap extension penalty [4]
    inputBinding:
      position: 1
      prefix: '-E'
  - id: '#R'
    type:
      - 'null'
      - int
    description: |
      INT    stop searching when there are >INT equally best hits [30]
    inputBinding:
      position: 1
      prefix: '-R'
  - id: '#q'
    type:
      - 'null'
      - int
    description: |
      INT    quality threshold for read trimming down to 35bp [0]
    inputBinding:
      position: 1
      prefix: '-q'
  - id: '#f'
    type:
      - 'null'
      - File
    description: |
      FILE   file to write output to instead of stdout
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#B'
    type:
      - 'null'
      - int
    description: |
      INT    length of barcode
    inputBinding:
      position: 1
      prefix: '-B'
  - id: '#L'
    type:
      - 'null'
      - boolean
    description: "       log-scaled gap penalty for long deletions\n"
    inputBinding:
      position: 1
      prefix: '-L'
  - id: '#N'
    type:
      - 'null'
      - boolean
    description: "       non-iterative mode: search for all n-difference hits (slooow)\n"
    inputBinding:
      position: 1
      prefix: '-N'
  - id: '#I'
    type:
      - 'null'
      - boolean
    description: "       the input is in the Illumina 1.3+ FASTQ-like format\n"
    inputBinding:
      position: 1
      prefix: '-I'
  - id: '#b'
    type:
      - 'null'
      - boolean
    description: "       the input read file is in the BAM format\n"
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#0'
    type:
      - 'null'
      - boolean
    description: "       use single-end reads only (effective with -b)\n"
    inputBinding:
      position: 1
      prefix: '-0'
  - id: '#1'
    type:
      - 'null'
      - boolean
    description: "       use the 1st read in a pair (effective with -b)\n"
    inputBinding:
      position: 1
      prefix: '-1'
  - id: '#2'
    type:
      - 'null'
      - boolean
    description: "       use the 2nd read in a pair (effective with -b)\n"
    inputBinding:
      position: 1
      prefix: '-2'
  - id: '#Y'
    type:
      - 'null'
      - boolean
    description: "       filter Casava-filtered sequences\n"
    inputBinding:
      position: 1
      prefix: '-Y'
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
  - aln
description: |-
  Usage:   bwa aln [options] <prefix> <in.fq>

  Options: -n NUM    max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
           -o INT    maximum number or fraction of gap opens [1]
           -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]
           -i INT    do not put an indel within INT bp towards the ends [5]
           -d INT    maximum occurrences for extending a long deletion [10]
           -l INT    seed length [32]
           -k INT    maximum differences in the seed [2]
           -m INT    maximum entries in the queue [2000000]
           -t INT    number of threads [1]
           -M INT    mismatch penalty [3]
           -O INT    gap open penalty [11]
           -E INT    gap extension penalty [4]
           -R INT    stop searching when there are >INT equally best hits [30]
           -q INT    quality threshold for read trimming down to 35bp [0]
           -f FILE   file to write output to instead of stdout
           -B INT    length of barcode
           -L        log-scaled gap penalty for long deletions
           -N        non-iterative mode: search for all n-difference hits (slooow)
           -I        the input is in the Illumina 1.3+ FASTQ-like format
           -b        the input read file is in the BAM format
           -0        use single-end reads only (effective with -b)
           -1        use the 1st read in a pair (effective with -b)
           -2        use the 2nd read in a pair (effective with -b)
           -Y        filter Casava-filtered sequences


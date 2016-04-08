#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: samtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#filebam'
    type: File
    description: file.bam
    inputBinding:
      position: 2
  - id: '#c'
    type:
      - 'null'
      - int
    description: >
      --coverage <int>,<int>,<int>    Coverage distribution min,max,step
      [1,1000,1]
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#d'
    type:
      - 'null'
      - boolean
    description: >
      --remove-dups                   Exclude from statistics reads marked as
      duplicates
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#f'
    type:
      - 'null'
      - int
    description: >
      --required-flag  <str|int>      Required flag, 0 for unset. See also
      `samtools flags` [0]
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#F'
    type:
      - 'null'
      - float
    description: >
      --filtering-flag <str|int>      Filtering flag, 0 for unset. See also
      `samtools flags` [0]

      --GC-depth <float>              the size of GC-depth bins (decreasing bin
      size increases memory requirement) [2e4]
    inputBinding:
      position: 1
      prefix: '-F'
  - id: '#h'
    type:
      - 'null'
      - boolean
    description: |
      --help                          This help message
    inputBinding:
      position: 1
      prefix: '-h'
  - id: '#i'
    type:
      - 'null'
      - int
    description: |
      --insert-size <int>             Maximum insert size [8000]
    inputBinding:
      position: 1
      prefix: '-i'
  - id: '#I'
    type:
      - 'null'
      - string
    description: >
      --id <string>                   Include only listed read group or sample
      name
    inputBinding:
      position: 1
      prefix: '-I'
  - id: '#l'
    type:
      - 'null'
      - int
    description: >
      --read-length <int>             Include in the statistics only reads with
      the given read length []
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#m'
    type:
      - 'null'
      - float
    description: >
      --most-inserts <float>          Report only the main part of inserts
      [0.99]
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#q'
    type:
      - 'null'
      - int
    description: |
      --trim-quality <int>            The BWA trimming parameter [0]
    inputBinding:
      position: 1
      prefix: '-q'
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: >
      --ref-seq <file>                Reference sequence (required for GC-depth
      and mismatches-per-cycle calculation).
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#t'
    type:
      - 'null'
      - boolean
    description: >
      --target-regions <file>         Do stats in these regions only.
      Tab-delimited file chr,from,to, 1-based, inclusive.
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      --sam                           Input is SAM (usually auto-detected now).
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#x'
    type:
      - 'null'
      - boolean
    description: >
      --sparse                        Suppress outputting IS rows where there are
      no insertions.
    inputBinding:
      position: 1
      prefix: '-x'
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
  - samtools
  - stats
description: |-
  About: The program collects statistics from BAM files. The output can be visualized using plot-bamstats.
  Usage: samtools stats [OPTIONS] file.bam
         samtools stats [OPTIONS] file.bam chr:from-to
  Options:
      -c, --coverage <int>,<int>,<int>    Coverage distribution min,max,step [1,1000,1]
      -d, --remove-dups                   Exclude from statistics reads marked as duplicates
      -f, --required-flag  <str|int>      Required flag, 0 for unset. See also `samtools flags` [0]
      -F, --filtering-flag <str|int>      Filtering flag, 0 for unset. See also `samtools flags` [0]
          --GC-depth <float>              the size of GC-depth bins (decreasing bin size increases memory requirement) [2e4]
      -h, --help                          This help message
      -i, --insert-size <int>             Maximum insert size [8000]
      -I, --id <string>                   Include only listed read group or sample name
      -l, --read-length <int>             Include in the statistics only reads with the given read length []
      -m, --most-inserts <float>          Report only the main part of inserts [0.99]
      -q, --trim-quality <int>            The BWA trimming parameter [0]
      -r, --ref-seq <file>                Reference sequence (required for GC-depth and mismatches-per-cycle calculation).
      -t, --target-regions <file>         Do stats in these regions only. Tab-delimited file chr,from,to, 1-based, inclusive.
      -s, --sam                           Input is SAM (usually auto-detected now).
      -x, --sparse                        Suppress outputting IS rows where there are no insertions.


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
  - id: '#[inbam]'
    type:
      - 'null'
      - File
    description: '[in.bam]'
    inputBinding:
      position: 3
  - id: '#]'
    type: boolean
    description: '...]'
    inputBinding:
      position: 2
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      INT     Set compression level, from 0 (uncompressed) to 9 (best)
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#m'
    type:
      - 'null'
      - int
    description: |
      INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#n'
    type:
      - 'null'
      - boolean
    description: "        Sort by read name\n"
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#o'
    type:
      - 'null'
      - File
    description: |
      FILE    Write final output to FILE rather than standard output
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#O'
    type:
      - 'null'
      - boolean
    description: |
      FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')   (either -O or
    inputBinding:
      position: 1
      prefix: '-O'
  - id: '#T'
    type:
      - 'null'
      - int
    description: |
      PREFIX  Write temporary files to PREFIX.nnnn.bam       -T is required)
      -@ INT     Set number of sorting and compression threads [1]
      Legacy usage: samtools sort [options...] <in.bam> <out.prefix>
    inputBinding:
      position: 1
      prefix: '-T'
  - id: '#f'
    type:
      - 'null'
      - boolean
    description: "        Use <out.prefix> as full final filename rather than prefix\n"
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#o'
    type:
      - 'null'
      - boolean
    description: "        Write final output to stdout rather than <out.prefix>.bam\n-l,m,n,@   Similar to corresponding options above\n"
    inputBinding:
      position: 1
      prefix: '-o'
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
  - sort
description: |-
  Usage: samtools sort [options...] [in.bam]
  Options:
    -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
    -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
    -n         Sort by read name
    -o FILE    Write final output to FILE rather than standard output
    -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')   (either -O or
    -T PREFIX  Write temporary files to PREFIX.nnnn.bam       -T is required)
    -@ INT     Set number of sorting and compression threads [1]

  Legacy usage: samtools sort [options...] <in.bam> <out.prefix>
  Options:
    -f         Use <out.prefix> as full final filename rather than prefix
    -o         Write final output to stdout rather than <out.prefix>.bam
    -l,m,n,@   Similar to corresponding options above


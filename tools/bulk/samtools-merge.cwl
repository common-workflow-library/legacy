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
  - id: '#in2bam'
    type: File
    description: '<in2.bam>'
    inputBinding:
      position: 7
  - id: '#in1bam'
    type: File
    description: '<in1.bam>'
    inputBinding:
      position: 6
  - id: '#outbam'
    type: File
    description: '<out.bam>'
    inputBinding:
      position: 5
  - id: '#n'
    type:
      - 'null'
      - boolean
    description: |
      sort by read names
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: "      attach RG tag (inferred from file names)\n"
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: "      uncompressed BAM output\n"
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#f'
    type:
      - 'null'
      - boolean
    description: "      overwrite the output BAM if exist\n"
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#1'
    type:
      - 'null'
      - boolean
    description: "      compress level 1\n"
    inputBinding:
      position: 1
      prefix: '-1'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      INT   compression level, from 0 to 9 [-1]
      -@ INT   number of BAM compression threads [0]
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#R'
    type:
      - 'null'
      - string
    description: |
      STR   merge file in the specified region STR [all]
    inputBinding:
      position: 1
      prefix: '-R'
  - id: '#h'
    type:
      - 'null'
      - File
    description: |
      FILE  copy the header in FILE to <out.bam> [in1.bam]
    inputBinding:
      position: 1
      prefix: '-h'
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: "      combine RG tags with colliding IDs rather than amending them\n"
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#p'
    type:
      - 'null'
      - boolean
    description: "      combine PG tags with colliding IDs rather than amending them\n"
    inputBinding:
      position: 1
      prefix: '-p'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      VALUE override random seed
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#b'
    type:
      - 'null'
      - File
    description: |
      FILE  list of input BAM filenames, one per line [null]
    inputBinding:
      position: 1
      prefix: '-b'
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
  - merge
description: |-
  Usage:   samtools merge [-nurlf] [-h inh.sam] [-b <bamlist.fofn>] <out.bam> <in1.bam> <in2.bam> [<in3.bam> ... <inN.bam>]

  Options: -n       sort by read names
           -r       attach RG tag (inferred from file names)
           -u       uncompressed BAM output
           -f       overwrite the output BAM if exist
           -1       compress level 1
           -l INT   compression level, from 0 to 9 [-1]
           -@ INT   number of BAM compression threads [0]
           -R STR   merge file in the specified region STR [all]
           -h FILE  copy the header in FILE to <out.bam> [in1.bam]
           -c       combine RG tags with colliding IDs rather than amending them
           -p       combine PG tags with colliding IDs rather than amending them
           -s VALUE override random seed
           -b FILE  list of input BAM filenames, one per line [null]


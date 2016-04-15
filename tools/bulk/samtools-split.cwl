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
  - id: '#]'
    type: boolean
    description: ']'
    inputBinding:
      position: 3
  - id: '#f'
    type:
      - 'null'
      - string
    description: |
      STRING       output filename format string ["%*_%#.bam"]
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: |
      FILE1        put reads with no RG tag or an unrecognised RG tag in FILE1
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: |
      FILE1:FILE2  ...and override the header with FILE2
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#v'
    type:
      - 'null'
      - string
    description: |
                   verbose output
      Format string expansions:
      %%     %
      %*     basename
      %#     @RG index
      %!     @RG ID
    inputBinding:
      position: 1
      prefix: '-v'
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
  - split
description: |-
  Usage: samtools split [-u <unaccounted.bam>[:<unaccounted_header.sam>]]
                        [-f <format_string>] [-v] <merged.bam>
  Options:
    -f STRING       output filename format string ["%*_%#.bam"]
    -u FILE1        put reads with no RG tag or an unrecognised RG tag in FILE1
    -u FILE1:FILE2  ...and override the header with FILE2
    -v              verbose output

  Format string expansions:
    %%     %
    %*     basename
    %#     @RG index
    %!     @RG ID


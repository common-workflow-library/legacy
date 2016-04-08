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
  - id: '#inbaminsamincram'
    type: File
    description: '<in.bam>|<in.sam>|<in.cram>'
    inputBinding:
      position: 2
  - id: '#b'
    type:
      - 'null'
      - boolean
    description: |
      output BAM
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#C'
    type:
      - 'null'
      - boolean
    description: "      output CRAM (requires -T)\n"
    inputBinding:
      position: 1
      prefix: '-C'
  - id: '#1'
    type:
      - 'null'
      - boolean
    description: "      use fast BAM compression (implies -b)\n"
    inputBinding:
      position: 1
      prefix: '-1'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: "      uncompressed BAM output (implies -b)\n"
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#h'
    type:
      - 'null'
      - boolean
    description: "      include header in SAM output\n"
    inputBinding:
      position: 1
      prefix: '-h'
  - id: '#H'
    type:
      - 'null'
      - boolean
    description: "      print SAM header only (no alignments)\n"
    inputBinding:
      position: 1
      prefix: '-H'
  - id: '#c'
    type:
      - 'null'
      - boolean
    description: "      print only the count of matching records\n"
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#o'
    type:
      - 'null'
      - File
    description: |
      FILE  output file name [stdout]
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#U'
    type:
      - 'null'
      - File
    description: |
      FILE  output reads not selected by filters to FILE [null]
    inputBinding:
      position: 1
      prefix: '-U'
  - id: '#t'
    type:
      - 'null'
      - boolean
    description: |
      FILE  FILE listing reference names and lengths (see long help) [null]
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#T'
    type:
      - 'null'
      - File
    description: |
      FILE  reference sequence FASTA FILE [null]
    inputBinding:
      position: 1
      prefix: '-T'
  - id: '#L'
    type:
      - 'null'
      - File
    description: |
      FILE  only include reads overlapping this BED FILE [null]
    inputBinding:
      position: 1
      prefix: '-L'
  - id: '#r'
    type:
      - 'null'
      - string
    description: |
      STR   only include reads in read group STR [null]
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#R'
    type:
      - 'null'
      - File
    description: |
      FILE  only include reads with read group listed in FILE [null]
    inputBinding:
      position: 1
      prefix: '-R'
  - id: '#q'
    type:
      - 'null'
      - int
    description: |
      INT   only include reads with mapping quality >= INT [0]
    inputBinding:
      position: 1
      prefix: '-q'
  - id: '#l'
    type:
      - 'null'
      - string
    description: |
      STR   only include reads in library STR [null]
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#m'
    type:
      - 'null'
      - int
    description: |
      INT   only include reads with number of CIGAR operations
      consuming query sequence >= INT [0]
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#f'
    type:
      - 'null'
      - int
    description: |
      INT   only include reads with all bits set in INT set in FLAG [0]
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#F'
    type:
      - 'null'
      - int
    description: |
      INT   only include reads with none of the bits set in INT
      set in FLAG [0]
    inputBinding:
      position: 1
      prefix: '-F'
  - id: '#x'
    type:
      - 'null'
      - string
    description: |
      STR   read tag to strip (repeatable) [null]
    inputBinding:
      position: 1
      prefix: '-x'
  - id: '#B'
    type:
      - 'null'
      - boolean
    description: "      collapse the backward CIGAR operation\n"
    inputBinding:
      position: 1
      prefix: '-B'
  - id: '#s'
    type:
      - 'null'
      - float
    description: |
      FLOAT integer part sets seed of random number generator [0];
      rest sets fraction of templates to subsample [no subsampling]
      -@ INT   number of BAM compression threads [0]
      -?       print long help, including note about region specification
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: "      ignored (input format is auto-detected)\n"
    inputBinding:
      position: 1
      prefix: '-S'
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
  - view
description: |-
  Usage:   samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]

  Options: -b       output BAM
           -C       output CRAM (requires -T)
           -1       use fast BAM compression (implies -b)
           -u       uncompressed BAM output (implies -b)
           -h       include header in SAM output
           -H       print SAM header only (no alignments)
           -c       print only the count of matching records
           -o FILE  output file name [stdout]
           -U FILE  output reads not selected by filters to FILE [null]
           -t FILE  FILE listing reference names and lengths (see long help) [null]
           -T FILE  reference sequence FASTA FILE [null]
           -L FILE  only include reads overlapping this BED FILE [null]
           -r STR   only include reads in read group STR [null]
           -R FILE  only include reads with read group listed in FILE [null]
           -q INT   only include reads with mapping quality >= INT [0]
           -l STR   only include reads in library STR [null]
           -m INT   only include reads with number of CIGAR operations
                    consuming query sequence >= INT [0]
           -f INT   only include reads with all bits set in INT set in FLAG [0]
           -F INT   only include reads with none of the bits set in INT
                    set in FLAG [0]
           -x STR   read tag to strip (repeatable) [null]
           -B       collapse the backward CIGAR operation
           -s FLOAT integer part sets seed of random number generator [0];
                    rest sets fraction of templates to subsample [no subsampling]
           -@ INT   number of BAM compression threads [0]
           -?       print long help, including note about region specification
           -S       ignored (input format is auto-detected)


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
  - id: '#[in2fq]'
    type:
      - 'null'
      - File
    description: '[in2.fq]'
    inputBinding:
      position: 4
  - id: '#in1fq'
    type: File
    description: '<in1.fq>'
    inputBinding:
      position: 3
  - id: '#idxbase'
    type: boolean
    description: '<idxbase>'
    inputBinding:
      position: 2
  - id: '#t'
    type:
      - 'null'
      - int
    description: |
      INT        number of threads [1]
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#k'
    type:
      - 'null'
      - int
    description: |
      INT        minimum seed length [19]
    inputBinding:
      position: 1
      prefix: '-k'
  - id: '#w'
    type:
      - 'null'
      - int
    description: |
      INT        band width for banded alignment [100]
    inputBinding:
      position: 1
      prefix: '-w'
  - id: '#d'
    type:
      - 'null'
      - int
    description: |
      INT        off-diagonal X-dropoff [100]
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#r'
    type:
      - 'null'
      - float
    description: >
      FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT
      [1.5]
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#y'
    type:
      - 'null'
      - int
    description: |
      INT        seed occurrence for the 3rd round seeding [20]
    inputBinding:
      position: 1
      prefix: '-y'
  - id: '#c'
    type:
      - 'null'
      - int
    description: |
      INT        skip seeds with more than INT occurrences [500]
    inputBinding:
      position: 1
      prefix: '-c'
  - id: '#D'
    type:
      - 'null'
      - float
    description: >
      FLOAT      drop chains shorter than FLOAT fraction of the longest
      overlapping chain [0.50]
    inputBinding:
      position: 1
      prefix: '-D'
  - id: '#W'
    type:
      - 'null'
      - int
    description: |
      INT        discard a chain if seeded bases shorter than INT [0]
    inputBinding:
      position: 1
      prefix: '-W'
  - id: '#m'
    type:
      - 'null'
      - int
    description: |
      INT        perform at most INT rounds of mate rescues for each read [50]
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: "           skip mate rescue\n"
    inputBinding:
      position: 1
      prefix: '-S'
  - id: '#P'
    type:
      - 'null'
      - boolean
    description: "           skip pairing; mate rescue performed unless -S also in use\n"
    inputBinding:
      position: 1
      prefix: '-P'
  - id: '#e'
    type:
      - 'null'
      - boolean
    description: "           discard full-length exact matches\nScoring options:\n"
    inputBinding:
      position: 1
      prefix: '-e'
  - id: '#A'
    type:
      - 'null'
      - int
    description: >
      INT        score for a sequence match, which scales options -TdBOELU unless
      overridden [1]
    inputBinding:
      position: 1
      prefix: '-A'
  - id: '#B'
    type:
      - 'null'
      - int
    description: |
      INT        penalty for a mismatch [4]
    inputBinding:
      position: 1
      prefix: '-B'
  - id: '#O'
    type:
      - 'null'
      - int
    description: |
      INT[,INT]  gap open penalties for deletions and insertions [6,6]
    inputBinding:
      position: 1
      prefix: '-O'
  - id: '#E'
    type:
      - 'null'
      - int
    description: >
      INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k'
      [1,1]
    inputBinding:
      position: 1
      prefix: '-E'
  - id: '#L'
    type:
      - 'null'
      - int
    description: |
      INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
    inputBinding:
      position: 1
      prefix: '-L'
  - id: '#U'
    type:
      - 'null'
      - int
    description: |
      INT        penalty for an unpaired read pair [17]
    inputBinding:
      position: 1
      prefix: '-U'
  - id: '#x'
    type:
      - 'null'
      - string
    description: |
      STR        read type. Setting -x changes multiple parameters unless overriden [null]
      pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
      ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
      intractg: -B9 -O16 -L5  (intra-species contigs to ref)
      Input/output options:
    inputBinding:
      position: 1
      prefix: '-x'
  - id: '#p'
    type:
      - 'null'
      - boolean
    description: "           smart pairing (ignoring in2.fq)\n"
    inputBinding:
      position: 1
      prefix: '-p'
  - id: '#R'
    type:
      - 'null'
      - string
    description: "STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n"
    inputBinding:
      position: 1
      prefix: '-R'
  - id: '#H'
    type:
      - 'null'
      - string
    description: >
      STR/FILE   insert STR to header if it starts with @; or insert lines in
      FILE [null]
    inputBinding:
      position: 1
      prefix: '-H'
  - id: '#j'
    type:
      - 'null'
      - boolean
    description: "           treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n"
    inputBinding:
      position: 1
      prefix: '-j'
  - id: '#v'
    type:
      - 'null'
      - int
    description: |
      INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
    inputBinding:
      position: 1
      prefix: '-v'
  - id: '#T'
    type:
      - 'null'
      - int
    description: |
      INT        minimum score to output [30]
    inputBinding:
      position: 1
      prefix: '-T'
  - id: '#h'
    type:
      - 'null'
      - int
    description: >
      INT[,INT]  if there are <INT hits with score >80% of the max score, output
      all in XA [5,200]
    inputBinding:
      position: 1
      prefix: '-h'
  - id: '#a'
    type:
      - 'null'
      - boolean
    description: "           output all alignments for SE or unpaired PE\n"
    inputBinding:
      position: 1
      prefix: '-a'
  - id: '#C'
    type:
      - 'null'
      - boolean
    description: "           append FASTA/FASTQ comment to SAM output\n"
    inputBinding:
      position: 1
      prefix: '-C'
  - id: '#V'
    type:
      - 'null'
      - boolean
    description: "           output the reference FASTA header in the XR tag\n"
    inputBinding:
      position: 1
      prefix: '-V'
  - id: '#Y'
    type:
      - 'null'
      - boolean
    description: "           use soft clipping for supplementary alignments\n"
    inputBinding:
      position: 1
      prefix: '-Y'
  - id: '#M'
    type:
      - 'null'
      - boolean
    description: "           mark shorter split hits as secondary\n"
    inputBinding:
      position: 1
      prefix: '-M'
  - id: '#I'
    type:
      - 'null'
      - float
    description: |
      FLOAT[,FLOAT[,INT[,INT]]]
      specify the mean, standard deviation (10% of the mean if absent), max
      (4 sigma from the mean if absent) and min of the insert size distribution.
      FR orientation only. [inferred]
      Note: Please read the man page for detailed description of the command line and options.
    inputBinding:
      position: 1
      prefix: '-I'
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
  - mem
description: |
  Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]

  Algorithm options:

         -t INT        number of threads [1]
         -k INT        minimum seed length [19]
         -w INT        band width for banded alignment [100]
         -d INT        off-diagonal X-dropoff [100]
         -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
         -y INT        seed occurrence for the 3rd round seeding [20]
         -c INT        skip seeds with more than INT occurrences [500]
         -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
         -W INT        discard a chain if seeded bases shorter than INT [0]
         -m INT        perform at most INT rounds of mate rescues for each read [50]
         -S            skip mate rescue
         -P            skip pairing; mate rescue performed unless -S also in use
         -e            discard full-length exact matches

  Scoring options:

         -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
         -B INT        penalty for a mismatch [4]
         -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
         -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
         -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
         -U INT        penalty for an unpaired read pair [17]

         -x STR        read type. Setting -x changes multiple parameters unless overriden [null]
                       pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
                       ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
                       intractg: -B9 -O16 -L5  (intra-species contigs to ref)

  Input/output options:

         -p            smart pairing (ignoring in2.fq)
         -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
         -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
         -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)

         -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
         -T INT        minimum score to output [30]
         -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
         -a            output all alignments for SE or unpaired PE
         -C            append FASTA/FASTQ comment to SAM output
         -V            output the reference FASTA header in the XR tag
         -Y            use soft clipping for supplementary alignments
         -M            mark shorter split hits as secondary

         -I FLOAT[,FLOAT[,INT[,INT]]]
                       specify the mean, standard deviation (10% of the mean if absent), max
                       (4 sigma from the mean if absent) and min of the insert size distribution.
                       FR orientation only. [inferred]

  Note: Please read the man page for detailed description of the command line and options.


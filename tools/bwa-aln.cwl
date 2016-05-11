#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

requirements:
  - $import: envvar-global.yml
  - class: InlineJavascriptRequirement

hints:
  - $import: bwa-docker.yml

inputs:
  - id: "prefix"
    type: File
    inputBinding:
      position: 4

  - id: "input"
    type: File
    inputBinding:
      position: 5

  - id: "output_filename"
    type: string
    inputBinding:
      position: 1
      prefix: "-f"

  - id: "threads"
    type: ["null",int]
    inputBinding:
      position: 1
      prefix: "-t"

  - id: "n"
    type: ["null",int,float]
    description: |
      max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
    inputBinding:
      position: 1
      prefix: "-n"

  - id: "o"
    type: ["null",int]
    description: |
      maximum number or fraction of gap opens [1]
    inputBinding:
      position: 1
      prefix: "-o"

  - id: "e"
    type: ["null",int]
    description: |
      maximum number of gap extensions, -1 for disabling long gaps [-1]
    inputBinding:
      position: 1
      prefix: "-e"

  - id: "i"
    type: ["null",int]
    description: |
      do not put an indel within INT bp towards the ends [5]
    inputBinding:
      position: 1
      prefix: "-i"

  - id: "d"
    type: ["null",int]
    description: |
      maximum occurrences for extending a long deletion [10]
    inputBinding:
      position: 1
      prefix: "-d"

  - id: "l"
    type: ["null",int]
    description: |
      seed length [32]
    inputBinding:
      position: 1
      prefix: "-l"

  - id: "k"
    type: ["null",int]
    description: |
      maximum differences in the seed [2]
    inputBinding:
      position: 1
      prefix: "-k"

  - id: "m"
    type: ["null",int]
    description: |
      maximum entries in the queue [2000000]
    inputBinding:
      position: 1
      prefix: "-m"

  - id: "I"
    type: ["null",boolean]
    description: |
      the input is in the Illumina 1.3+ FASTQ-like format
    inputBinding:
      position: 1
      prefix: "-I"

outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

baseCommand:
  - bwa
  - aln

description: |
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
           -c        input sequences are in the color space
           -L        log-scaled gap penalty for long deletions
           -N        non-iterative mode: search for all n-difference hits (slooow)
           -I        the input is in the Illumina 1.3+ FASTQ-like format
           -b        the input read file is in the BAM format
           -0        use single-end reads only (effective with -b)
           -1        use the 1st read in a pair (effective with -b)
           -2        use the 2nd read in a pair (effective with -b)

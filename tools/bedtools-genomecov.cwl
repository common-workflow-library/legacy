#!/usr/bin/env cwl-runner

class: CommandLineTool

description: |
  Tool:    bedtools genomecov (aka genomeCoverageBed)
  Sources: https://github.com/arq5x/bedtools2
  Summary: Compute the coverage of a feature file among a genome.
  Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>

requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: bedtools-docker.cwl
  - import: bedtools-genomecov-types.cwl

inputs:
  - id: "#input"
    type: File
    description: |
      The input file can be in BAM format
          (Note: BAM _must_ be sorted by position)
      or <bed/gff/vcf>
    inputBinding:
      position: 1
      secondaryFiles:
        - engine: node-engine.cwl
          script: |
           {
            if ((/.*\.bam$/i).test($job['input'].path))
               return {"path": $job['input'].path+".bai", "class": "File"};
            return [];
           }
      valueFrom:
        engine: node-engine.cwl
        script: |
          {
            var prefix = ((/.*\.bam$/i).test($job['input'].path))?'-ibam':'-i';
            return [prefix,$job['input'].path];
          }

  - id: "#genomeFile"
    type: File
    description:
      Input genome file.
    inputBinding:
      position: 2
      prefix: "-g"

  - id: "#dept"
    type: ["null","#depts"]

  - id: "#scale"
    type: ["null",float ]
    description: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
    inputBinding:
      position: 4
      prefix: -scale

  - id: "#dz"
    type: boolean
    description: |
      Report the depth at each genome position (with zero-based coordinates).
      Reports only non-zero positions.
      Default behavior is to report a histogram.
    default: false
    inputBinding:
      position: 4
      prefix: "-dz"

  - id: "#split"
    type: boolean
    description: |
      reat "split" BAM or BED12 entries as distinct BED intervals.
      when computing coverage.
      For BAM files, this uses the CIGAR "N" and "D" operations
      to infer the blocks for computing coverage.
      For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
      fields (i.e., columns 10,11,12).
    default: false
    inputBinding:
      position: 4
      prefix: "-split"

  - id: "#strand"
    type: ["null", string]
    description: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6).
      - (STRING): can be + or -
    inputBinding:
      position: 4
      prefix: "-strand"

  - id: "#max"
    type: ["null",int]
    description: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
    inputBinding:
      position: 4
      prefix: "-max"

  - id: "#m5"
    type: boolean
    description: |
      Calculate coverage of 5" positions (instead of entire interval).
    default: false
    inputBinding:
      position: 4
      prefix: "-5"

  - id: "#m3"
    type: boolean
    description: |
      Calculate coverage of 3" positions (instead of entire interval).
    default: false
    inputBinding:
      position: 4
      prefix: "-3"

  - id: "#genomecoverageout"
    type: string

outputs:
  - id: "#genomecoverage"
    type: File
    description: "The file containing the genome coverage"
    outputBinding:
      glob: 
        engine: cwl:JsonPointer
        script: /job/genomecoverageout
stdout: 
  engine: cwl:JsonPointer
  script: /job/genomecoverageout

baseCommand: ["bedtools", "genomecov"]

#!/usr/bin/env cwl-runner
class: CommandLineTool

requirements:
  - import: node-engine-local.cwl

#hints:
#  - class: DockerRequirement
#    dockerPull: 
#    dockerImageId: 

description: |
  Tool:    bedtools genomecov (aka genomeCoverageBed)
  Version: v2.24.0
  Summary: Compute the coverage of a feature file among a genome.

  Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>

inputs:
  - id: "#input"
    type: File
    description: |
      The input file is in BAM format.
      Note: BAM _must_ be sorted by position
      Or <bed/gff/vcf>
    inputBinding:
      secondaryFiles:
        - engine: node-engine-local.cwl
          script: |
           {
            if ((/.*\.bam$/i).test($job.inputs['input'].path))
               return {"path": $job.inputs['input'].path+".bai", "class": "File"};
            return [];
           }
      
  - id: "#genomeFile"
    type: File
    description:
      Input genome file.
    inputBinding:
      prefix: "-g"

  - id: "#scale"
    type: ["null",float ]
    description: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
    inputBinding:
      position: 1
      prefix: -scale

  - id: "#d"
    type: boolean
    description: |
      Report the depth at each genome position (with one-based coordinates).
      Default behavior is to report a histogram.
    default: false
    inputBinding:
      position: 1
      prefix: "-d"

  - id: "#dz"
    type: boolean
    description: |
      Report the depth at each genome position (with zero-based coordinates).
      Reports only non-zero positions.
      Default behavior is to report a histogram.
    default: false
    inputBinding:
      position: 1
      prefix: "-dz"

  - id: "#bg"
    type: boolean
    description: |
      Report depth in BedGraph format. For details, see
      genome.ucsc.edu/goldenPath/help/bedgraph.html
    default: false
    inputBinding:
      position: 1
      prefix: "-bg"
      valueFrom:
        engine: node-engine-local.cwl
        script: |
          {
           if ($job.inputs['scale'] == null) {
              return null;
            } else {
              return true;
            }
          }

  - id: "#bga"
    type: boolean
    description: |
      Report depth in BedGraph format, as above (-bg).
      However with this option, regions with zero 
      coverage are also reported. This allows one to
      quickly extract all regions of a genome with 0 
      coverage by applying: "grep -w 0$" to the output.
    default: false
    inputBinding:
      position: 1
      prefix: "-bga"

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
      position: 1
      prefix: "-split"

  - id: "#strand"
    type: ["null", string]
    description: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6). 
      - (STRING): can be + or -
    inputBinding:
      position: 1
      prefix: "-strand"

  - id: "#max"
    type: ["null",int]
    description: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
    inputBinding:
      position: 1
      prefix: "-max"

  - id: "#m5"
    type: boolean
    description: |
      Calculate coverage of 5" positions (instead of entire interval).
    default: false
    inputBinding:
      position: 1
      prefix: "-5"

  - id: "#m3"
    type: boolean
    description: |
      Calculate coverage of 3" positions (instead of entire interval).
    default: false
    inputBinding:
      position: 1
      prefix: "-3"

outputs:
  - id: "#genecoverage"
    type: File
    description: "The file containing the genome coverage"
    outputBinding:
      glob: genomecov.out
stdout: genomecov.out

baseCommand: ["bedtools", "genomecov"]

arguments:
  - valueFrom:
      engine: node-engine-local.cwl
      script: |
        {
          var param="-i";
          if ((/.*\.bam$/i).test($job.inputs['input'].path)) param="-ibam";
          return [param, $job.inputs['input'].path];
        }

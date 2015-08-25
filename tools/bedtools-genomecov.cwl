#!/usr/bin/env cwl-runner
#bedtools genomecov -bg -split -scale 0.03967752645856 -ibam "$1.bam" -g /wardrobe/indices/STAR/mm10/chrNameLength.txt

class: CommandLineTool

description: "Invoke 'bedtools genomecov' "

inputs:
  - id: "#input"
    type: File
    description: |
      The input file is in BAM format.
      Note: BAM _must_ be sorted by position
    inputBinding:
      prefix: "-ibam"
      secondaryFiles:
        - ".bai"

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

  - id: "#bg"
    type: boolean
    description: |
      Report depth in BedGraph format. For details, see
      genome.ucsc.edu/goldenPath/help/bedgraph.html
    default: false
    inputBinding:
      position: 1
      prefix: "-bg"
      #valueFrom:
      #  engine: "cwl:JsonPointer"
      #  script: |
      #    {
      #      if ($job['scale']) {
      #        return true;
      #      } else {
      #        return null;
      #      }
      #    }

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

outputs:
  - id: "#genecoverage"
    type: File
    description: "The file containing the genome coverage"
    outputBinding:
      glob: genomecov.out
stdout: genomecov.out

baseCommand: ["/usr/local/bin/bedtools", "genomecov"]

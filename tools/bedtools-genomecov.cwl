#!/usr/bin/env cwl-runner
#bedtools genomecov -bg -split -scale 0.03967752645856 -ibam "$1.bam" -g /wardrobe/indices/STAR/mm10/chrNameLength.txt

class: CommandLineTool

description: "Invoke 'bedtools genomecov' "

inputs:
  - id: "#input"
    type: File
    description:
      Input bam file.
    inputBinding:
      prefix: "-ibam"

  - id: "#genomeFile"
    type: File
    description:
      Input genome file.
    inputBinding:
      prefix: "-g"

outputs:
  - id: "#genecoverage"
    type: File
    description: "The file containing the genome coverage"
    outputBinding:
      glob: genomecoverage.out
stdout: genomecoverage.out

baseCommand: ["/usr/local/bin/bedtools", "genomecov"]

arguments:
  - "-bg"
  - "-split"

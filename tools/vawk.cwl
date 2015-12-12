#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3.dev2"

class: CommandLineTool

requirements:
 - "@import": envvar-global.cwl

description: |
 usage: vawk [-h] [-v VAR] [-c INFO_COL] [--header] [--debug] cmd [vcf]
 positional arguments:
   cmd                   vawk command syntax is exactly the same as awk syntax with
        a few additional features. The INFO field can be split using
        the I$ prefix and the SAMPLE field can be split using
        the S$ prefix. For example, I$AF prints the allele frequency of
        each variant and S$NA12878 prints the entire SAMPLE field for the
        NA12878 individual for each variant. S$* returns all samples.
        The SAMPLE field can be further split based on the keys in the
        FORMAT field of the VCF (column 9). For example, S$NA12877$GT
        returns the genotype of the NA12878 individual.
        ex: '{ if (I$AF>0.5) print $1,$2,$3,I$AN,S$NA12878,S$NA12877$GT }'
   vcf                   VCF file (default: stdin)
 optional arguments:
   -h, --help            show this help message and exit
   -v VAR, --var VAR     declare an external variable (e.g.: SIZE=10000)
   -c INFO_COL, --col INFO_COL
    column of the INFO field [8]
   --header              print VCF header
   --debug               debugging level verbosity

inputs:
  - id: "#cmd"
    type: string
    description: |
        vawk command syntax is exactly the same as awk syntax with a few
        additional features. The INFO field can be split using the I$ prefix
        and the SAMPLE field can be split using the S$ prefix. For example,
        I$AF prints the allele frequency of each variant and S$NA12878 prints
        the entire SAMPLE field for the NA12878 individual for each variant.
        S$* returns all samples. The SAMPLE field can be further split based on
        the keys in the FORMAT field of the VCF (column 9). For example,
        S$NA12877$GT returns the genotype of the NA12878 individual.
        ex: '{ if (I$AF>0.5) print $1,$2,$3,I$AN,S$NA12878,S$NA12877$GT }'
    inputBinding:
      position: 1
    streamable: true

  - id: "#input"
    type: File
    description: |
      VCF file
    inputBinding:
      position: 2

stdout:
   "output.vcf"

outputs:
  - id: "#processed"
    type: File
    description: "The resulting VCF file"
    streamable: true
    outputBinding:
      glob: "output.vcf"
      

baseCommand: ["vawk"]

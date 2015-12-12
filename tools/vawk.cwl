#!/usr/bin/env cwl-runner

"@context":
  "foaf": "http://xmlns.com/foaf/0.1/"
  "doap": "http://usefulinc.com/ns/doap"
  "adms": "http://purl.org/adms/"
  "admssw": "http://purl.org/adms/sw/"

adms:Asset:
  admssw:SoftwareProject:
    doap:name: "vawk"
    doap:description: >
      An awk-like VCF parser
    doap:homepage: "https://github.com/cc2qe/vawk"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/cc2qe/vawk.git"
    doap:release:
      - doap:revision: "0.0.1"
    doap:license: "None"
    doap:category: "commandline tool"
    doap:programming-language: "Python"
    doap:developer:
      - foaf:Person:
        foaf:name: "Colby Chiang"
        foaf:Organization: "Washington University"
  adms:AssetDistribution:
    doap:name: "vawk.cwl"
    doap:description: "Developed for CWL consortium http://commonwl.org/"
    doap:specification: "http://common-workflow-language.github.io/draft-3/"
    doap:release: "cwl:draft-3.dev2"
    doap:homepage: "http://commonwl.org/"
    doap:location: "https://github.com/common-workflow-language/workflows/blob/master/tools/vawk.cwl"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/common-workflow-language/workflows"
    doap:maintainer:
      foaf:Person:
        foaf:openid: "http://orcid.org/0000-0002-2961-9670"
        foaf:name: "Michael R. Crusoe"
        foaf:mbox: "mailto:crusoe@ucdavis.edu"
        foaf:organization: "University of California, Davis"

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

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

requirements:
  - "@import": envvar-global.cwl

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
  - id: "#sorted"
    type: File
    description: "The resulting VCF file"
    streamable: true
    outputBinding: "output.vcf"
      

baseCommand: ["vawk"]

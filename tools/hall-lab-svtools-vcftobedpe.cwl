#!/usr/bin/env cwl-runner

$namespaces:
  "cwl": "https://w3id.org/cwl/cwl#"
  "foaf": "http://xmlns.com/foaf/0.1/"
  "doap": "http://usefulinc.com/ns/doap"
  "adms": "http://purl.org/adms/"
  "admssw": "http://purl.org/adms/sw/"

$schemas:
  - https://joinup.ec.europa.eu/svn/adms_foss/adms_sw_v1.00/adms_sw_v1.00.rdf
  - http://xmlns.com/foaf/spec/20140114.rdf
  - http://usefulinc.com/ns/doap#
  - http://www.w3.org/ns/adms

adms:Asset:
  admssw:SoftwareProject:
    doap:name: "svtools"
    doap:description: >
      Comprehensive utilities to explore structural variations in genomes.
    doap:homepage: "https://github.com/hall-lab/svtools"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/hall-lab/svtools"
    doap:release:
      - doap:revision: "0.0.1-44188e60c44c4"
    doap:license: "no license"
    doap:category: "commandline tool"
    doap:programming-language: "Python"
  adms:AssetDistribution:
    doap:name: "samtools-index.cwl"
    doap:description: "Developed for CWL consortium http://commonwl.org/"
    doap:specification: "http://common-workflow-language.github.io/draft-3/"
    doap:release: "cwl:draft-3.dev2"
    doap:homepage: "http://commonwl.org/"
    doap:location: "https://github.com/common-workflow-language/workflows/blob/master/tools/hall-lab-svtools.cwl"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/common-workflow-language/workflows"

cwlVersion: "cwl:draft-3.dev2"

class: CommandLineTool

description: |
  Usage: vcftobedbpe -i <in.vcf> -o [out.bedpe]

requirements:
  - "@import": envvar-global.cwl

inputs:
  - id: "#input"
    type: File
    description: |
      "Input vcf file."
    inputBinding:
      prefix: "-i"


outputs:
  - id: "#output"
    type: File
    description: "The bedpe file"
    outputBinding:
      prefix: "-o"

baseCommand: ["samtools"]

arguments:
  - valueFrom: $(inputs.bai?'-b':inputs.csi?'-c':[])
    position: 1
  - valueFrom: $(new_ext())
    position: 3

#!/usr/bin/env cwl-runner

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  adms: http://www.w3.org/ns/adms#
  dcat: http://www.w3.org/ns/dcat#

$schemas:
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf

cwlVersion: "cwl:draft-3"

class: CommandLineTool

adms:includedAsset:
  doap:name: "sort"
  doap:description: >
    sort - sort lines of text files
  doap:release:
  - class: doap:Version
    doap:revision: "5.93"
  doap:license: "GPL"
  doap:category: "commandline tool"
  doap:programming-language: "C"
  doap:developer:
  - class: foaf:Person
    foaf:name: "Mike Haertel"
  - class: foaf:Person
    foaf:name: "Paul Eggert"
  doap:mailing-list:
    foaf:mbox: "bug-coreutils@gnu.org"

description: |
  linux-sort.cwl is developed for CWL consortium

doap:name: "linux-sort.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/linux-sort.cwl"

doap:maintainer:
- class: foaf:Organization
  foaf:name: "Barski Lab, Cincinnati Children's Hospital Medical Center"
  foaf:member:
  - class: foaf:Person
    id: "http://orcid.org/0000-0001-9102-5681"
    foaf:openid: "http://orcid.org/0000-0001-9102-5681"
    foaf:name: "Andrey Kartashov"
    foaf:mbox: "mailto:Andrey.Kartashov@cchmc.org"

requirements:
  - $import: envvar-global.yml
  - $import: linux-sort-docker.yml
  - class: InlineJavascriptRequirement

inputs:
  - id: "#input"
    type:
      type: array
      items: File
    inputBinding:
      position: 4

  - id: "#output"
    type: string

  - id: "#key"
    type: 
      type: array
      items: string
      inputBinding:
        prefix: "-k"
    inputBinding:
      position: 1
    description: |
      -k, --key=POS1[,POS2]
      start a key at POS1, end it at POS2 (origin 1)

outputs:
  - id: "#sorted"
    type: File
    description: "The sorted file"
    outputBinding:
      glob: $(inputs.output)

stdout: $(inputs.output)

baseCommand: ["sort"]



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

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

adms:includedAsset:
  $include: alea-ontology.yaml

description: |
  alea-phaseVCF.cwl is developed for CWL consortium

dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/alea-phaseVCF.cwl"
doap:repository:
- class: doap:GitRepository
  doap:location: "https://github.com/common-workflow-language/workflows"
doap:homepage: "http://commonwl.org/"
doap:license: "Apache2"

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
  - $import: envvar-global.cwl
  - $import: alea-docker.cwl
  - class: InlineJavascriptRequirement

inputs:
  - id: "#hapsDir"
    type: File
    description: |
      path to the directory containing the .haps files
    inputBinding:
      position: 2

  - id: "#unphased"
    type: File
    description: |
      path to the vcf file containing unphased SNPs and Indels
    inputBinding:
      position: 3

  - id: "#outputPrefix"
    type: string
    description: |
      output file prefix including the path but not the extension
    inputBinding:
      position: 3

outputs:
  - id: "#phasevcf"
    type: File
    description: "Creates the file outputPrefix.vcf.gz"
    outputBinding:
      glob: $(inputs.outputPrefix+".vcf.gz")

baseCommand: ["alea", "phaseVCF"]

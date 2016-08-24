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
  doap:name: "picard"
  doap:description: >
    A set of Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats.
    Picard is implemented using the HTSJDK Java library HTSJDK, supporting accessing of common file formats,
    such as SAM and VCF, used for high-throughput sequencing data.
    http://broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex
  doap:homepage: "http://broadinstitute.github.io/picard/"
  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/broadinstitute/picard.git"
  doap:release:
  - class: doap:Version
    doap:revision: "1.141"
  doap:license: "MIT, Apache2"
  doap:category: "commandline tool"
  doap:programming-language: "JAVA"
  doap:developer:
  - class: foaf:Organization
    foaf:name: "Broad Institute"

description: |
  picard-BuildBamIndex.cwl is developed for CWL consortium
    Generates a BAM index (.bai) file.

doap:name: "picard-BuildBamIndex.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/picard-BuildBamIndex.cwl"

dct:creator:
- class: foaf:Organization
  foaf:name: "UNIVERSITY OF MELBOURNE"
  foaf:member:
  - class: foaf:Person
    id: "farahk@student.unimelb.edu.au"
    foaf:mbox: "mailto:farahk@student.unimelb.edu.au"
  - class: foaf:Person
    id: "skanwal@student.unimelb.edu.au"
    foaf:mbox: "mailto:skanwal@student.unimelb.edu.au"

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
- $import: picard-docker.yml
- class: InlineJavascriptRequirement

inputs:

- id: "java_arg"
  type: string
  default: "-Xmx4g"
  inputBinding:
    position: 1

- id: "INPUT"
  type: File
  description: >
   INPUT String A BAM file or URL to process. Must be sorted in coordinate order.
  inputBinding:
    position: 4
    separate: false
    prefix: "INPUT="

outputs:
  - id: "index"
    type: File
    outputBinding: 
      glob: $(inputs.OUTPUT)

baseCommand: ["java"]

arguments:
- valueFrom: "/usr/local/bin/picard.jar"
  position: 2
  prefix: "-jar"
- valueFrom: "BuildBamIndex"
  position: 3
- valueFrom: $(inputs.INPUT.path.split('/').slice(-1)[0]+".bai")
  position: 5
  separate: false
  prefix: "OUTPUT="
#  description: >
#    OUTPUT File The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.
#    If INPUT is a URL and OUTPUT is unspecified, defaults to a file in the current directory.

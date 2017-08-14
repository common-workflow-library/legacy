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

cwlVersion: v1.0
class: CommandLineTool

adms:includedAsset:
  doap:name: picard
  doap:description: 'A set of Java command line tools for manipulating high-throughput
    sequencing data (HTS) data and formats. Picard is implemented using the HTSJDK
    Java library HTSJDK, supporting accessing of common file formats, such as SAM
    and VCF, used for high-throughput sequencing data. http://broadinstitute.github.io/picard/command-line-overview.html#SortSam

    '
  doap:homepage: http://broadinstitute.github.io/picard/
  doap:repository:
  - class: doap:GitRepository
    doap:location: https://github.com/broadinstitute/picard.git
  doap:release:
  - class: doap:Version
    doap:revision: '1.141'
  doap:license: MIT, Apache2
  doap:category: commandline tool
  doap:programming-language: JAVA
  doap:developer:
  - class: foaf:Organization
    foaf:name: Broad Institute
doap:name: picard-SortSam.cwl
dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/picard-SortSam.cwl
dct:creator:
- class: foaf:Organization
  foaf:name: THE UNIVERSITY OF MELBOURNE
  foaf:member:
  - class: foaf:Person
    id: farahk@student.unimelb.edu.au
    foaf:name: Farah Zaib Khan
    foaf:mbox: mailto:farahk@student.unimelb.edu.au
  - class: foaf:Person
    id: skanwal@student.unimelb.edu.au
    foaf:name: Sehrish Kanwal
    foaf:mbox: mailto:skanwal@student.unimelb.edu.au
requirements:
- $import: picard-docker.yml
- class: InlineJavascriptRequirement

inputs:
  inputFileName_sortSam:
    type: File
    inputBinding:
      position: 4
      prefix: INPUT=
    doc: The BAM or SAM file to sort. Required
  outputFileName_sortSam:
    type: string
    inputBinding:
      position: 5
      prefix: OUTPUT=
    doc: The sorted BAM or SAM output file. Required
  createIndex:
    type: string?
    default: 'true'
    inputBinding:
      position: 8
      prefix: CREATE_INDEX=
    doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
      Default value false. This option can be set to 'null' to clear the default value.
      Possible values {true, false}
  SO-coordinate:
    type: string
    default: coordinate
    inputBinding:
      position: 6
      prefix: SORT_ORDER=
    doc: Sort order of output file Required. Possible values {unsorted, queryname,
      coordinate, duplicate}
  java_arg:
    type: string
    default: -Xmx2g
    inputBinding:
      position: 1

  tmpdir:
    type: string
    inputBinding:
      position: 7
      prefix: TMP_DIR=
    doc: Default value null. This option may be specified 0 or more times.
outputs:
  sortSam_output:
    type: File
    outputBinding:
      glob: $(inputs.outputFileName_sortSam)

arguments:
- valueFrom: /usr/local/bin/picard.jar
  position: 2
  prefix: -jar
- valueFrom: SortSam
  position: 3

baseCommand: [java]
doc: |
  picard-SortSam.cwl is developed for CWL consortium
    Generates a sorted file


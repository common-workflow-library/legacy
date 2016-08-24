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
    http://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary
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
  picard-CreateSequenceDictionary.cwl is developed for CWL consortium
  Read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary.

doap:name: "picard-CreateSequenceDictionary.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/picard-CreateSequenceDictionary.cwl"

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
- id: "reference"
  type: File
  description: |
    Input reference fasta or fasta.gz
  inputBinding:
    prefix: "REFERENCE="
    separate: false
    position: 4

- id: "output_filename"
  type: string
  description: |
    Output SAM or BAM file containing only the sequence dictionary
  inputBinding:
    prefix: "OUTPUT="
    separate: false
    position: 4

- id: "GENOME_ASSEMBLY"
  type: ["null",string]
  description: |
    Put into AS field of sequence dictionary entry if supplied
  inputBinding:
    prefix: "GENOME_ASSEMBLY="
    separate: false
    position: 4

- id: "URI"
  type: ["null",string]
  description: |
    Put into UR field of sequence dictionary entry.
    If not supplied, input reference file is used
  inputBinding:
    prefix: "URI="
    separate: false
    position: 4

- id: "SPECIES"
  type: ["null",string]
  description: |
    Put into SP field of sequence dictionary entry
  inputBinding:
    prefix: "SPECIES="
    separate: false
    position: 4

- id: "TRUNCATE_NAMES_AT_WHITESPACE"
  type: ["null",boolean]
  description: |
    Make sequence name the first word from the > line in the fasta file.  By default the
    entire contents of the > line is used, excluding leading and trailing whitespace.
    Default value: true. This option can be set to 'null' to clear the default value.
    Possible values: {true, false}
  inputBinding:
    prefix: "TRUNCATE_NAMES_AT_WHITESPACE="
    separate: false
    position: 4

- id: "NUM_SEQUENCES"
  type: ["null",int]
  description: |
    Stop after writing this many sequences.  For testing.
    Default value: 2147483647.
  inputBinding:
    prefix: "NUM_SEQUENCES="
    separate: false
    position: 4

outputs:
- id: "#output"
  type: File
  description: >
    Output SAM or BAM file containing only the sequence dictionary Required.
  outputBinding:
    glob: $(inputs.output_filename)

baseCommand: ["java"]

arguments:
- valueFrom: "-Xmx4g"
  position: 1
- valueFrom: "/usr/local/bin/picard.jar"
  position: 2
  prefix: "-jar"
- valueFrom: "CreateSequenceDictionary"
  position: 3


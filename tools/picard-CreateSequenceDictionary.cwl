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
    and VCF, used for high-throughput sequencing data. http://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary

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
doap:name: picard-CreateSequenceDictionary.cwl
dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/picard-CreateSequenceDictionary.cwl
doap:maintainer:
- class: foaf:Organization
  foaf:name: Barski Lab, Cincinnati Children's Hospital Medical Center
  foaf:member:
  - class: foaf:Person
    id: http://orcid.org/0000-0001-9102-5681
    foaf:openid: http://orcid.org/0000-0001-9102-5681
    foaf:name: Andrey Kartashov
    foaf:mbox: mailto:Andrey.Kartashov@cchmc.org
requirements:
- $import: picard-docker.yml
- class: InlineJavascriptRequirement

inputs:
  reference:
    type: File
    inputBinding:
      prefix: REFERENCE=
      separate: false
      position: 4

    doc: |
      Input reference fasta or fasta.gz
  URI:
    type: string?
    inputBinding:
      prefix: URI=
      separate: false
      position: 4

    doc: |
      Put into UR field of sequence dictionary entry.
      If not supplied, input reference file is used
  output_filename:
    type: string
    inputBinding:
      prefix: OUTPUT=
      separate: false
      position: 4

    doc: |
      Output SAM or BAM file containing only the sequence dictionary
  NUM_SEQUENCES:
    type: int?
    inputBinding:
      prefix: NUM_SEQUENCES=
      separate: false
      position: 4

    doc: |
      Stop after writing this many sequences.  For testing.
      Default value: 2147483647.
  TRUNCATE_NAMES_AT_WHITESPACE:
    type: boolean?
    inputBinding:
      prefix: TRUNCATE_NAMES_AT_WHITESPACE=
      separate: false
      position: 4

    doc: |
      Make sequence name the first word from the > line in the fasta file.  By default the
      entire contents of the > line is used, excluding leading and trailing whitespace.
      Default value: true. This option can be set to 'null' to clear the default value.
      Possible values: {true, false}
  GENOME_ASSEMBLY:
    type: string?
    inputBinding:
      prefix: GENOME_ASSEMBLY=
      separate: false
      position: 4

    doc: |
      Put into AS field of sequence dictionary entry if supplied
  SPECIES:
    type: string?
    inputBinding:
      prefix: SPECIES=
      separate: false
      position: 4

    doc: |
      Put into SP field of sequence dictionary entry
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

    doc: 'Output SAM or BAM file containing only the sequence dictionary Required.

      '
baseCommand: [java]
arguments:
- valueFrom: -Xmx4g
  position: 1
- valueFrom: /usr/local/bin/picard.jar
  position: 2
  prefix: -jar
- valueFrom: CreateSequenceDictionary
  position: 3

doc: |
  picard-CreateSequenceDictionary.cwl is developed for CWL consortium
  Read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary.


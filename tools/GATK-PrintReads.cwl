#!/usr/bin/env cwl-runner

#$namespaces:
#  dct: http://purl.org/dc/terms/
#  foaf: http://xmlns.com/foaf/0.1/
#  doap: http://usefulinc.com/ns/doap#
#  adms: http://www.w3.org/ns/adms#
#  dcat: http://www.w3.org/ns/dcat#
#
#$schemas:
#- http://dublincore.org/2012/06/14/dcterms.rdf
#- http://xmlns.com/foaf/spec/20140114.rdf
#- http://usefulinc.com/ns/doap#
#- http://www.w3.org/ns/adms#
#- http://www.w3.org/ns/dcat.rdf

cwlVersion: v1.0
class: CommandLineTool

#adms:includedAsset:
#  doap:name: GATK
#  doap:description: 'The Genome Analysis Toolkit or GATK is a software package for
#    analysis of high-throughput sequencing data, developed by the Data Science and
#    Data Engineering group at the Broad Institute.  The toolkit offers a wide variety
#    of tools, with a primary focus on variant discovery and genotyping as well as
#    strong emphasis on data quality assurance. Its robust architecture, powerful processing
#    engine and high-performance computing features make it capable of taking on projects
#    of any size. http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles
#
#    '
#  doap:homepage: https://www.broadinstitute.org/gatk/
#  doap:repository:
#  - class: doap:GitRepository
#    doap:location: https://github.com/broadgsa/gatk.git
#  doap:release:
#  - class: doap:Version
#    doap:revision: '3.4'
#  doap:license: mixed licensing model
#  doap:category: commandline tool
#  doap:programming-language: JAVA
#  doap:developer:
#  - class: foaf:Organization
#    foaf:name: Broad Institute
#doap:name: GATK-PrintReads.cwl
#dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-PrintReads.cwl
#dct:creator:
#- class: foaf:Organization
#  foaf:name: THE UNIVERSITY OF MELBOURNE
#  foaf:member:
#  - class: foaf:Person
#    id: farahk@student.unimelb.edu.au
#    foaf:name: Farah Zaib Khan
#    foaf:mbox: mailto:farahk@student.unimelb.edu.au
#  - class: foaf:Person
#    id: skanwal@student.unimelb.edu.au
#    foaf:name: Sehrish Kanwal
#    foaf:mbox: mailto:skanwal@student.unimelb.edu.au
#doap:maintainer:
#- class: foaf:Organization
#  foaf:name: THE UNIVERSITY OF MELBOURNE
#  foaf:member:
#  - class: foaf:Person
#    id: farahk@student.unimelb.edu.au
#    foaf:name: Farah Zaib Khan
#    foaf:mbox: mailto:farahk@student.unimelb.edu.au
#  - class: foaf:Person
#    id: skanwal@student.unimelb.edu.au
#    foaf:name: Sehrish Kanwal
#    foaf:mbox: mailto:skanwal@student.unimelb.edu.au
requirements:
- $import: envvar-global.yml
- $import: envvar-global.yml
- $import: GATK-docker.yml

inputs:
  sample_file:
    type: File[]?
    inputBinding:
      position: 11



  reference:
    type: File
    inputBinding:
      position: 5
      prefix: -R
    secondaryFiles:
    - .64.amb
    - .64.ann
    - .64.bwt
    - .64.pac
    - .64.sa
    - .fai
    - ^.dict
  platform:
    type: string?
    inputBinding:
      position: 13
      prefix: --platform
    doc: Exclude all reads with this platform from the output
  input_baseRecalibrator:
    type: File
    inputBinding:
      position: 7
      prefix: -BQSR
    doc: the recalibration table produced by BaseRecalibration
  number:
    type: string?
    inputBinding:
      position: 13
      prefix: --number
    doc: Exclude all reads with this platform from the output
  simplify:
    type: boolean?
    inputBinding:
      position: 9
      prefix: --simplify
    doc: Erase all extra attributes in the read but keep the read group information
  readGroup:
    type: string?
    inputBinding:
      position: 12
      prefix: --readGroup
    doc: Exclude all reads with this read group from the output
  sample_name:
    type: string[]?
    inputBinding:
      position: 10

    doc: Sample name to be included in the analysis. Can be specified multiple times.
  inputBam_printReads:
    type: File
    inputBinding:
      position: 6
      prefix: -I
    secondaryFiles:
    - ^.bai
    doc: bam file produced after indelRealigner
  outputfile_printReads:
    type: string?
    inputBinding:
      position: 8
      prefix: -o
    doc: name of the output file from indelRealigner
  java_arg:
    type: string
    default: -Xmx4g
    inputBinding:
      position: 1

  threads:
    type: int
    default: 4
    inputBinding:
      prefix: -nct
      position: 14

outputs:
  output_printReads:
    type: File
    outputBinding:
      glob: $(inputs.outputfile_printReads)


arguments:
- valueFrom: ./test/test-files
  position: 2
  separate: false
  prefix: -Djava.io.tmpdir=
- valueFrom: /usr/local/bin/GenomeAnalysisTK.jar
  position: 3
  prefix: -jar
- valueFrom: PrintReads
  position: 4
  prefix: -T
baseCommand: [java]
doc: |
  GATK-RealignTargetCreator.cwl is developed for CWL consortium
  Prints all reads that have a mapping quality above zero
    Usage: java -Xmx4g -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta -I input1.bam -I input2.bam -o output.bam --read_filter MappingQualityZero

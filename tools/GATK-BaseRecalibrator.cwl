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
  doap:name: GATK
  doap:description: 'The Genome Analysis Toolkit or GATK is a software package for
    analysis of high-throughput sequencing data, developed by the Data Science and
    Data Engineering group at the Broad Institute.  The toolkit offers a wide variety
    of tools, with a primary focus on variant discovery and genotyping as well as
    strong emphasis on data quality assurance. Its robust architecture, powerful processing
    engine and high-performance computing features make it capable of taking on projects
    of any size. http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles

    '
  doap:homepage: https://www.broadinstitute.org/gatk/
  doap:repository:
  - class: doap:GitRepository
    doap:location: https://github.com/broadgsa/gatk.git
  doap:release:
  - class: doap:Version
    doap:revision: '3.4'
  doap:license: mixed licensing model
  doap:category: commandline tool
  doap:programming-language: JAVA
  doap:developer:
  - class: foaf:Organization
    foaf:name: Broad Institute
doap:name: GATK-BaseRecalibrator.cwl
dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-BaseRecalibrator.cwl
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
doap:maintainer:
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
- $import: envvar-global.yml
- $import: GATK-docker.yml

inputs:
  deletions_default_quality:
    type: int?
    inputBinding:
      position: 24
      prefix: --deletions_default_quality
    doc: default quality for the base deletions covariate
  reference:
    type: File
    inputBinding:
      position: 5
      prefix: -R
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .pac
    - .rbwt
    - .rpac
    - .rsa
    - .sa
    - .fai
    - ^.dict
  binary_tag_name:
    type: string?
    inputBinding:
      position: 27
      prefix: --binary_tag_name
    doc: the binary tag covariate name if using it
  no_standard_covs:
    type: boolean?
    inputBinding:
      position: 15
      prefix: --no_standard_covs
    doc: Do not use the standard set of covariates, but rather just the ones listed
      using the -cov argument
  solid_nocall_strategy:
    type: string?
    inputBinding:
      position: 11
      prefix: --solid_nocall_strategy
    doc: Defines the behavior of the recalibrator when it encounters no calls in the
      color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ
  low_quality_tail:
    type: int?
    inputBinding:
      position: 20
      prefix: --low_quality_tail
    doc: minimum quality for the bases in the tail of the reads to be considered
  java_arg:
    type: string
    default: -Xmx4g
    inputBinding:
      position: 1

  out:
    type: File?
    inputBinding:
      position: 14
      prefix: --out
    doc: The output recalibration table file to create
  quantizing_levels:
    type: boolean?
    inputBinding:
      position: 13
      prefix: --quantizing_levels
    doc: Sort the rows in the tables of reports. Whether GATK report tables should
      have rows in sorted order, starting from leftmost column
  bqsrBAQGapOpenPenalty:
    type: double?
    inputBinding:
      position: 26
      prefix: --bqsrBAQGapOpenPenalty
    doc: BQSR BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is perhaps
      better for whole genome call sets
  mismatches_context_size:
    type: int?
    inputBinding:
      position: 17
      prefix: --mismatches_context_size
    doc: Size of the k-mer context to be used for base mismatches
  maximum_cycle_value:
    type: int?
    inputBinding:
      position: 18
      prefix: --maximum_cycle_value
    doc: The maximum cycle value permitted for the Cycle covariate
  run_without_dbsnp_potentially_ruining_quality:
    type: boolean?
    inputBinding:
      position: 12
      prefix: --run_without_dbsnp_potentially_ruining_quality
    doc: If specified, allows the recalibrator to be used without a dbsnp rod. Very
      unsafe and for expert users only.
  inputBam_BaseRecalibrator:
    type: File
    inputBinding:
      position: 6
      prefix: -I
    secondaryFiles:
    - ^.bai
    doc: bam file produced after indelRealigner
  known:
    type: File[]?
    inputBinding:
      position: 7

  outputfile_BaseRecalibrator:
    type: string
    inputBinding:
      position: 8
      prefix: -o
    doc: name of the output file from baseRecalibrator
  lowMemoryMode:
    type: boolean?
    inputBinding:
      position: 19
      prefix: --lowMemoryMode
    doc: Reduce memory usage in multi-threaded code at the expense of threading efficiency
  solid_recal_mode:
    type: string?
    inputBinding:
      position: 10
      prefix: --solid_recal_mode
    doc: How should we recalibrate solid bases in which the reference was inserted?
      Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS
  insertions_default_quality:
    type: int?
    inputBinding:
      position: 22
      prefix: --insertions_default_quality
    doc: default quality for the base insertions covariate
  sort_by_all_columns:
    type: boolean?
    inputBinding:
      position: 9
      prefix: --sort_by_all_columns
    doc: Sort the rows in the tables of reports. Whether GATK report tables should
      have rows in sorted order, starting from leftmost column
  list:
    type: boolean?
    inputBinding:
      position: 21
      prefix: --list
    doc: List the available covariates and exit
  indels_context_size:
    type: int?
    inputBinding:
      position: 23
      prefix: --indels_context_size
    doc: Size of the k-mer context to be used for base insertions and deletions
  mismatches_default_quality:
    type: int?
    inputBinding:
      position: 16
      prefix: --mismatches_default_quality
    doc: default quality for the base mismatches covariate
  covariate:
    type: string?
    inputBinding:
      position: 25
      prefix: --covariate
    doc: One or more covariates to be used in the recalibration. Can be specified
      multiple times
outputs:
  output_baseRecalibrator:
    type: File
    outputBinding:
      glob: $(inputs.outputfile_BaseRecalibrator)

arguments:
- valueFrom: ./test/test-files
  position: 2
  separate: false
  prefix: -Djava.io.tmpdir=
- valueFrom: /usr/local/bin/GenomeAnalysisTK.jar
  position: 3
  prefix: -jar
- valueFrom: BaseRecalibrator
  position: 4
  prefix: -T
baseCommand: [java]
doc: |
  GATK-BaseRecalibrator.cwl is developed for CWL consortium
  It generate base recalibration table to compensate for systematic errors in basecalling confidences
    Usage: java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fasta -I my_reads.bam -knownSites latest_dbsnp.vcf -o recal_data.table.


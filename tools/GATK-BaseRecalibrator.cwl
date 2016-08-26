#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool


hints:
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
    - .64.amb
    - .64.ann
    - .64.bwt
    - .64.pac
    - .64.sa
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
    type:
      - "null"
      - type: array
        items: File
        inputBinding:
          prefix: --knownSites
    inputBinding:
      position: 28
    doc: 'Any number of VCF files representing known SNPs and/or indels. Could be
      e.g. dbSNP and/or official 1000 Genomes indel calls. SNPs in these files will
      be ignored unless the --mismatchFraction argument is used. optional parameter.'

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
    type:
      type: array
      items: string
      inputBinding:
        prefix: --covariate
    inputBinding:
      position: 25
    doc: One or more covariates to be used in the recalibration. Can be specified
      multiple times

    threads:
      type: int
      default: 4
      inputBinding:
        prefix: -nct
        position: 26

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

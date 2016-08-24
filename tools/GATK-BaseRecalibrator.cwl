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
  doap:name: "GATK"
  doap:description: >
    The Genome Analysis Toolkit or GATK is a software package for analysis of
    high-throughput sequencing data, developed by the Data Science and Data
    Engineering group at the Broad Institute. 
    The toolkit offers a wide variety of tools, with a primary focus on variant
    discovery and genotyping as well as strong emphasis on data quality
    assurance.
    Its robust architecture, powerful processing engine and high-performance
    computing features make it capable of taking on projects of any size.
    http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles
  doap:homepage: "https://www.broadinstitute.org/gatk/"
  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/broadgsa/gatk.git"
  doap:release:
  - class: doap:Version
    doap:revision: "3.4"
  doap:license: "mixed licensing model"
  doap:category: "commandline tool"
  doap:programming-language: "JAVA"
  doap:developer:
  - class: foaf:Organization
    foaf:name: "Broad Institute"

description: |
  GATK-BaseRecalibrator.cwl is developed for CWL consortium
  It generate base recalibration table to compensate for systematic errors in basecalling confidences
    Usage: java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fasta -I my_reads.bam -knownSites latest_dbsnp.vcf -o recal_data.table.
doap:name: "GATK-BaseRecalibrator.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-BaseRecalibrator.cwl"

dct:creator:
- class: foaf:Organization
  foaf:name: "THE UNIVERSITY OF MELBOURNE"
  foaf:member:
  - class: foaf:Person
    id: "farahk@student.unimelb.edu.au"
    foaf:name: "Farah Zaib Khan"
    foaf:mbox: "mailto:farahk@student.unimelb.edu.au"
    
  - class: foaf:Person
    id: "skanwal@student.unimelb.edu.au"
    foaf:name: "Sehrish Kanwal"
    foaf:mbox: "mailto:skanwal@student.unimelb.edu.au"

doap:maintainer:
- class: foaf:Organization
  foaf:name: "THE UNIVERSITY OF MELBOURNE"
  foaf:member:
  - class: foaf:Person
    id: "farahk@student.unimelb.edu.au"
    foaf:name: "Farah Zaib Khan"
    foaf:mbox: "mailto:farahk@student.unimelb.edu.au"
  - class: foaf:Person
    id: "skanwal@student.unimelb.edu.au"
    foaf:name: "Sehrish Kanwal"
    foaf:mbox: "mailto:skanwal@student.unimelb.edu.au"
    
requirements:
- $import: envvar-global.yml
- $import: GATK-docker.yml

inputs:
  - id: "#java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding: 
      position: 1
  
  - id: "#reference"
    type: File
    inputBinding:
      position: 5
      prefix: "-R"
    secondaryFiles:
      - ".amb"
      - ".ann"
      - ".bwt"
      - ".pac"
      - ".rbwt"
      - ".rpac"
      - ".rsa"
      - ".sa"
      - ".fai"
      - "^.dict"

  - id: "#inputBam_BaseRecalibrator"
    type: File
    description: bam file produced after indelRealigner
    inputBinding:
      position: 6 
      prefix: "-I"
    secondaryFiles:
      - "^.bai"

  - id: "#known"  
    type:
      - "null"
      - type: array
        items: File
        inputBinding: { prefix: "-knownSites" }
    inputBinding: 
      position: 7
      
  - id: "#outputfile_BaseRecalibrator"
    type: string
    description: name of the output file from baseRecalibrator
    inputBinding:
      position: 8
      prefix: "-o"

  - id: "#sort_by_all_columns"
    type: ["null", boolean]
    description: Sort the rows in the tables of reports. Whether GATK report tables should have rows in sorted order, starting from leftmost column
    inputBinding:
      position: 9
      prefix: "--sort_by_all_columns"

  - id: "#solid_recal_mode"
    type: ["null", string]
    description: How should we recalibrate solid bases in which the reference was inserted? Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS
    inputBinding:
      position: 10
      prefix: "--solid_recal_mode"

  - id: "#solid_nocall_strategy"
    type: ["null", string]
    description: Defines the behavior of the recalibrator when it encounters no calls in the color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ
    inputBinding:
      position: 11
      prefix: "--solid_nocall_strategy"


  - id: "#run_without_dbsnp_potentially_ruining_quality"
    type: ["null", boolean]
    description: If specified, allows the recalibrator to be used without a dbsnp rod. Very unsafe and for expert users only.
    inputBinding:
      position: 12
      prefix: "--run_without_dbsnp_potentially_ruining_quality"

  - id: "#quantizing_levels"
    type: ["null", boolean]
    description: Sort the rows in the tables of reports. Whether GATK report tables should have rows in sorted order, starting from leftmost column
    inputBinding:
      position: 13
      prefix: "--quantizing_levels"

  - id: "#out"
    type: ["null", File]
    description: The output recalibration table file to create
    inputBinding:
      position: 14
      prefix: "--out"

  - id: "#no_standard_covs"
    type: ["null", boolean]
    description: Do not use the standard set of covariates, but rather just the ones listed using the -cov argument
    inputBinding:
      position: 15
      prefix: "--no_standard_covs"

  - id: "#mismatches_default_quality"
    type: ["null", int]
    description: default quality for the base mismatches covariate
    inputBinding:
      position: 16
      prefix: "--mismatches_default_quality"
      
  - id: "#mismatches_context_size"
    type: ["null", int]
    description: Size of the k-mer context to be used for base mismatches
    inputBinding:
      position: 17
      prefix: "--mismatches_context_size"

  - id: "#maximum_cycle_value"
    type: ["null", int]
    description: The maximum cycle value permitted for the Cycle covariate
    inputBinding:
      position: 18
      prefix: "--maximum_cycle_value"

  - id: "#lowMemoryMode"
    type: ["null", boolean]
    description: Reduce memory usage in multi-threaded code at the expense of threading efficiency
    inputBinding:
      position: 19
      prefix: "--lowMemoryMode"
      
  - id: "#low_quality_tail"
    type: ["null", int]
    description: minimum quality for the bases in the tail of the reads to be considered
    inputBinding:
      position: 20
      prefix: "--low_quality_tail"
      
  - id: "#list"
    type: ["null", boolean]
    description: List the available covariates and exit
    inputBinding:
      position: 21
      prefix: "--list"

  - id: "#insertions_default_quality"
    type: ["null", int]
    description: default quality for the base insertions covariate
    inputBinding:
      position: 22
      prefix: "--insertions_default_quality"

  - id: "#indels_context_size"
    type: ["null", int]
    description: Size of the k-mer context to be used for base insertions and deletions
    inputBinding:
      position: 23
      prefix: "--indels_context_size"

  - id: "#deletions_default_quality"
    type: ["null", int]
    description: default quality for the base deletions covariate
    inputBinding:
      position: 24
      prefix: "--deletions_default_quality"

  - id: "#covariate"
    type: ["null", string]
    description: One or more covariates to be used in the recalibration. Can be specified multiple times
    inputBinding:
      position: 25
      prefix: "--covariate"

  - id: "#bqsrBAQGapOpenPenalty"
    type: ["null", double]
    description: BQSR BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is perhaps better for whole genome call sets
    inputBinding:
      position: 26
      prefix: "--bqsrBAQGapOpenPenalty"

  - id: "#binary_tag_name"
    type: ["null", string]
    description: the binary tag covariate name if using it
    inputBinding:
      position: 27
      prefix: "--binary_tag_name"

outputs:
  - id: "#output_baseRecalibrator"
    type: File
    outputBinding: 
      glob: $(inputs.outputfile_BaseRecalibrator)

arguments:

  - valueFrom: "./test/test-files"
    position: 2
    separate: false
    prefix: "-Djava.io.tmpdir="
    
  - valueFrom: "/usr/local/bin/GenomeAnalysisTK.jar"
    position: 3
    prefix: "-jar"

  - valueFrom: "BaseRecalibrator"
    position: 4
    prefix: "-T"

baseCommand: ["java"]


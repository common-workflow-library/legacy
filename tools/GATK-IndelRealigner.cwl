#!/usr/bin/env cwl-runner

#$namespaces:
##  dct: http://purl.org/dc/terms/
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
#
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
#doap:name: GATK-IndelRealigner.cwl
#dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-IndelRealigner.cwl
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
  inputBam_realign:
    type: File
    inputBinding:
      position: 6
      prefix: -I
    secondaryFiles:
    - ^.bai
    doc: bam file produced after markDups execution
  maxReadsForConsensuses:
    type: int?
    inputBinding:
      prefix: --maxReadsForConsensuses
    doc: Max reads used for finding the alternate consensuses (necessary to improve
      performance in deep coverage)
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
  LODThresholdForCleaning:
    type: double?
    inputBinding:
      prefix: --LODThresholdForCleaning
    doc: LOD threshold above which the cleaner will clean
  maxConsensuses:
    type: int?
    inputBinding:
      prefix: --maxConsensuses
    doc: Max alternate consensuses to try (necessary to improve performance in deep
      coverage)
  outputfile_indelRealigner:
    type: string
    inputBinding:
      position: 8
      prefix: -o
    doc: name of the output file from indelRealigner
  maxReadsInMemory:
    type: int?
    inputBinding:
      prefix: --maxReadsInMemory
    doc: max reads allowed to be kept in memory at a time by the SAMFileWriter
  maxIsizeForMovement:
    type: int?
    inputBinding:
      prefix: --maxIsizeForMovement
    doc: maximum insert size of read pairs that we attempt to realign. For expert
      users only!
  maxPositionalMoveAllowed:
    type: int?
    inputBinding:
      prefix: --maxPositionalMoveAllowed
    doc: Maximum positional move in basepairs that a read can be adjusted during realignment.
      For expert users only!
  bamout:
    type: File?
    inputBinding:
      prefix: --out
    doc: The realigned bam file. Optional parameter
  intervals:
    type: File
    inputBinding:
      position: 7
      prefix: -targetIntervals
    doc: list of intervals created by realignerTargetCreataor
  entropyThreshold:
    type: double?
    inputBinding:
      prefix: --entropyThreshold
    doc: Percentage of mismatches at a locus to be considered having high entropy
      (0.0 < entropy <= 1.0)
  maxReadsForRealignment:
    type: int?
    inputBinding:
      prefix: --maxReadsForRealignment
    doc: Max reads allowed at an interval for realignment
  known:
    type: File[]?
    doc: Any number of VCF files representing known SNPs and/or indels. Could be e.g.
      dbSNP and/or official 1000 Genomes indel calls. SNPs in these files will be
      ignored unless the --mismatchFraction argument is used. optional parameter.
  consensusDeterminationModel:
    type: string?
    inputBinding:
      prefix: --consensusDeterminationModel
    doc: Percentage of mismatches at a locus to be considered having high entropy
      (0.0 < entropy <= 1.0)
  noOriginalAlignmentTags:
    type: boolean?
    inputBinding:
      prefix: --noOriginalAlignmentTags
    doc: Dont output the original cigar or alignment start tags for each realigned
      read in the output bam
  java_arg:
    type: string
    default: -Xmx4g
    inputBinding:
      position: 1
# doesn't work, see http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
#  threads:
#    type: int
#    default: 4
#    inputBinding:
#      prefix: -nt
#      position: 4

  nWayOut:
    type: string?
    inputBinding:
      prefix: '--nWayOut '
    doc: Generate one output file for each input (-I) bam file (not compatible with
      -output). See the main page for more details.
outputs:
  output_indelRealigner:
    type: File
    outputBinding:
      glob: $(inputs.outputfile_indelRealigner)

arguments:
- valueFrom: ./test/test-files
  position: 2
  separate: false
  prefix: -Djava.io.tmpdir=
- valueFrom: /usr/local/bin/GenomeAnalysisTK.jar
  position: 3
  prefix: -jar
- valueFrom: IndelRealigner
  position: 4
  prefix: -T
baseCommand: [java]
doc: |
  GATK-RealignTargetCreator.cwl is developed for CWL consortium
  It performs local realignment of reads around indels.
    Usage: java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I input.bam --known indels.vcf -o forIndelRealigner.intervals.

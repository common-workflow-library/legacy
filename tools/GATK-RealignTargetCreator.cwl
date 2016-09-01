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
doap:name: GATK-RealignTargetCreator.cwl
dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-RealignTargetCreator.cwl
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

inputs: # position 0, for java args, 1 for the jar, 2 for the tool itself
  GATKJar:
    type: File
    inputBinding:
      position: 1
      prefix: "-jar"
  maxIntervalSize:
    type: int?
    inputBinding:
      prefix: --maxIntervalSize
      position: 2
    doc: 'maximum interval size; any intervals larger than this value will be dropped.
      optional paramter'
  inputBam_realign:
    type: File
    inputBinding:
      position: 2
      prefix: -I
    secondaryFiles:
    - ^.bai
    doc: 'bam file produced after mark-duplicates execution'
  outputfile_realignTarget:
    type: string
    inputBinding:
      prefix: -o
      position: 2
    doc: 'name of the output file from realignTargetCreator'
  reference:
    type: File
    inputBinding:
      position: 2
      prefix: -R
    secondaryFiles:
    - .64.amb
    - .64.ann
    - .64.bwt
    - .64.pac
    - .64.sa
    - .fai
    - ^.dict
    doc: 'human reference sequence along with the secondary files.'
  minReadsAtLocus:
    type: int?
    inputBinding:
      prefix: --minReadsAtLocus
      position: 2
    doc: 'minimum reads at a locus to enable using the entropy calculation'
  windowSize:
    type: int?
    inputBinding:
      prefix: --windowSize
      position: 2
    doc: 'window size for calculating entropy or SNP clusters'
  mismatchFraction:
    type: int?
    inputBinding:
      prefix: --mismatchFraction
      position: 2
    doc: 'fraction of base qualities needing to mismatch for a position to have high
      entropy'
  known:
    type:
      - "null" 
      - type: array
        items: File
        inputBinding:
          prefix: --known
    inputBinding:
      position: 2
    doc: 'Any number of VCF files representing known SNPs and/or indels. Could be
      e.g. dbSNP and/or official 1000 Genomes indel calls. SNPs in these files will
      be ignored unless the --mismatchFraction argument is used. optional parameter.'
  java_arg:
    type: string
    default: -Xmx4g
    inputBinding:
      position: 0

outputs:
  output_realignTarget:
    type: File
    outputBinding:
      glob: $(inputs.outputfile_realignTarget)

arguments:
- valueFrom: $(runtime.tmpdir)
  position: 0
  separate: false
  prefix: -Djava.io.tmpdir=
- valueFrom: RealignerTargetCreator
  position: 2
  prefix: -T
baseCommand: [java]
doc: |
  GATK-RealignTargetCreator.cwl is developed for CWL consortium
    It accepts 3 input files and produces a file containing list of target intervals to pass to the IndelRealigner.
    Usage: java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I input.bam --known indels.vcf -o forIndelRealigner.intervals.


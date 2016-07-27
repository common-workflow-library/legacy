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
    The Genome Analysis Toolkit or GATK is a software package for analysis of high-throughput sequencing data, developed by the Data Science and Data Engineering group at the Broad Institute. 
    The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance.
    Its robust architecture, powerful processing engine and high-performance computing features make it capable of taking on projects of any size.
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
  GATK-RealignTargetCreator.cwl is developed for CWL consortium
    It accepts 3 input files and produces a file containing list of target intervals to pass to the IndelRealigner.
    Usage: java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I input.bam --known indels.vcf -o forIndelRealigner.intervals.

doap:name: "GATK-RealignTargetCreator.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-RealignTargetCreator.cwl"

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
    description: >
      human reference sequence along with the secondary files.
    inputBinding:
      position: 5
      prefix: "-R"
    secondaryFiles:
      - ".amb"
      - ".ann"
      - ".bwt"
      - ".pac"
      - ".sa"
      - ".fai"
      - "^.dict"

  - id: "#inputBam_realign"
    type: File
    description: >
      bam file produced after mark-duplicates execution
    inputBinding:
      position: 6
      prefix: "-I"
    secondaryFiles:
      - "^.bai"

  - id: "#outputfile_realignTarget"
    type: string
    description: >
      name of the output file from realignTargetCreator
    inputBinding:
      prefix: "-o"
      position: 7

  - id: "#known"  
    type:
      - "null"
      - type: array
        items: File
        inputBinding: { prefix: "--known" }
    description: >
      Any number of VCF files representing known SNPs and/or indels. Could be e.g. dbSNP and/or official 1000 Genomes indel calls. SNPs in these files will be ignored unless the --mismatchFraction argument is used. optional parameter.
    inputBinding:
      position: 8

  - id: "#maxIntervalSize"  
    type: ["null", int]
    description: >
      maximum interval size; any intervals larger than this value will be dropped. optional paramter
    inputBinding:
      prefix: "--maxIntervalSize"

  - id: "#minReadsAtLocus"  
    type: ["null", int]
    description: >
      minimum reads at a locus to enable using the entropy calculation
    inputBinding:
      prefix: "--minReadsAtLocus"

  - id: "#mismatchFraction"  
    type: ["null", int]
    description: >
      fraction of base qualities needing to mismatch for a position to have high entropy
    inputBinding:
      prefix: "--mismatchFraction"
    
  - id: "#windowSize"  
    type: ["null", int]
    description: >
      window size for calculating entropy or SNP clusters
    inputBinding:
      prefix: "--windowSize"


outputs:
  - id: "#output_realignTarget"
    type: File
    outputBinding: 
      glob: $(inputs.outputfile_realignTarget)
        
arguments:
  - valueFrom: "./test/test-files"
    position: 2
    separate: false
    prefix: "-Djava.io.tmpdir="
    
  - valueFrom: "/usr/local/bin/GenomeAnalysisTK.jar"
    position: 3
    prefix: "-jar"

  - valueFrom: "RealignerTargetCreator"
    position: 4
    prefix: "-T"

baseCommand: ["java"]

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

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

adms:includedAsset:
  doap:name: "picard"
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

dct:isPartOf:
  doap:name: "CWL Workflows"
  doap:homepage: "http://commonwl.org/"
  doap:license: "Apache2"

  doap:implements:
  - class: doap:Specification
    doap:homepage: "http://common-workflow-language.github.io/draft-3/"

  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/common-workflow-language/workflows"

  dct:creator:
  - class: foaf:Organization
    foaf:name: "Curoverse"
    foaf:member:
    - class: foaf:Person
      id: "http://orcid.org/0000-0003-3566-7705"
      foaf:name: "Peter Amstutz"
      foaf:mbox: "mailto:peter.amstutz@curoverse.com"
  - class: foaf:Organization
    foaf:name: "Seven Bridges Genomics"
    foaf:member:
    - class: foaf:Person
      id: "mailto:nebojsa.tijanic@sbgenomics.com"
      foaf:name: "Nebojša Tijanić"
      foaf:mbox: "mailto:nebojsa.tijanic@sbgenomics.com"

  dct:contributor:
  - class: foaf:Organization
    foaf:name: "Seven Bridges Genomics"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Luka Stojanovic"
      foaf:mbox: "mailto:luka.stojanovic@sbgenomics.com"
  - class: foaf:Organization
    foaf:name: "Galaxy Project, Pennsylvania State University"
    foaf:member:
    - class: foaf:Person
      foaf:name: "John Chilton"
      foaf:mbox: "mailto:jmchilton@gmail.com"
  - class: foaf:Organization
    foaf:name: "University of California, Davis"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Michael R. Crusoe"
      foaf:mbox: "mailto:crusoe@ucdavis.edu"
  - class: foaf:Organization
    foaf:name: "Institut Pasteur"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Hervé Ménager"
      foaf:mbox: "mailto:herve.menager@gmail.com"
  - class: foaf:Organization
    foaf:name: "BioDatomics"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Maxim Mikheev"
      foaf:mbox: "mailto:mikhmv@biodatomics.com"
  - class: foaf:Organization
    foaf:name: "University of Manchester"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Stian Soiland-Reyes"
      foaf:mbox: "mailto:soiland-reyes@cs.manchester.ac.uk"

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
  foaf:name: "Barski Lab, Cincinnati Children's Hospital Medical Center"
  foaf:member:
  - class: foaf:Person
    id: "http://orcid.org/0000-0001-9102-5681"
    foaf:openid: "http://orcid.org/0000-0001-9102-5681"
    foaf:name: "Andrey Kartashov"
    foaf:mbox: "mailto:Andrey.Kartashov@cchmc.org"
    
requirements:
- $import: envvar-global.cwl
- $import: envvar-global.cwl
- $import: GATK-docker.cwl

inputs:

  - id: "#java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding: 
      position: 1
     
  - id: "#RealignerTarget"
    type: string
    description: >
      tool used for this step from GATK jar
    default: "RealignerTargetCreator"
    inputBinding: { position: 4, prefix: "-T" }

     
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

  - id: "#known"  
    type:
      type: array
      items: File
      inputBinding: { prefix: "--known" }
    description: >
      File set of known indels 
    inputBinding: 
      position: 7

  - id: "#outputfile_realignTarget"
    type: string
    description: >
      name of the output file from realignTargetCreator
    inputBinding:
      position: 8
      prefix: "-o"

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

baseCommand: ["java"]

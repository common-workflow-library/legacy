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
    http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles
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
  picard-MergeSamFiles.cwl is developed for CWL consortium

doap:name: "picard-MergeSamFiles.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/picard-MergeSamFiles.cwl"

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
- $import: envvar-global.yml
- $import: picard-docker.yml
- class: InlineJavascriptRequirement
  
inputs:
  - id: "java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding:
      position: 1

  - id: "#outputFileName_mergedSam"
    type: string
    description: >
      SAM or BAM file to write merged result to Required
    inputBinding:
      position: 4
      prefix: "OUTPUT="

  - id: "#inputFileName_mergedSam"
    type:
      type: array
      items: File
      inputBinding: { prefix: "INPUT=" } 
    description: >
      SAM or BAM input file Default value null. This option must be specified at least 1 times
    inputBinding: { position: 5}

  - id: "#readSorted"
    type: ["null", boolean]
    description: >
      If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    inputBinding:
      position: 6
      prefix: "ASSUME_SORTED="

  - id: "#mergeSequenceDictionaries"
    type: ["null", boolean]
    description: >
      Merge the sequence dictionaries Default value false. This option can be set to null to clear the default value. Possible values {true, false}
    inputBinding:
      position: 7
      prefix: "MERGE_SEQUENCE_DICTIONARIES="
      
  - id: "#useThreading"
    type: ["null", boolean]
    description: >
      Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    inputBinding:
      position: 8
      prefix: "USE_THREADING="

  - id: "#comment"
    type: ["null", string]
    description: >
      Comment(s) to include in the merged output files header. Default value null. This option may be specified 0 or more times
    inputBinding:
      position: 9
      prefix: "COMMENT="

  - id: "#tmpdir"
    type: string
    description: >
      Default value null. This option may be specified 0 or more times.
    inputBinding:
      position: 10
      prefix: "TMP_DIR="

outputs:
  - id: "#mergeSam_output"
    type: File
    outputBinding: 
      glob: $(inputs.outputFileName_mergedSam)
        
arguments: 

- valueFrom: "/usr/local/bin/picard.jar"
  position: 2
  prefix: "-jar"
  
- valueFrom: "MergeSamFiles"
  position: 3

  
baseCommand: ["java"]

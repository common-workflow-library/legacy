#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/
#The tool works with docker from broadinstitute/picard and works without -jar and java args making the command simple. 

class: CommandLineTool

description: |
  Reads FASTA or FASTA.GZ files containing reference sequences, and writes them as a SAM file containing a sequence dictionary
  Options: 
    REFERENCE File Input reference fasta or fasta.gz Required
    OUTPUT File Output SAM or BAM file containing only the sequence dictionary Required
    GENOME_ASSEMBLY String Put into AS field of sequence dictionary entry if supplied Default value null
    URI String Put into UR field of sequence dictionary entry. If not supplied, input reference file is used Default value null
    SPECIES String Put into SP field of sequence dictionary entry Default value null
    TRUNCATE_NAMES_AT_WHITESPACE  Boolean Make sequence name the first word from the > line in the fasta file. By default the entire contents of the > line is used, excluding leading and trailing whitespace. Default value true. This option can be set to 'null' to clear the default value. Possible values {true, false}
    NUM_SEQUENCES Integer Stop after writing this many sequences. For testing. Default value 2147483647. This option can be set to 'null' to clear the default value

requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: picard-docker.cwl

inputs:

  - id: "#createDictionary"
    type: string
    default: "CreateSequenceDictionary"
    inputBinding:
      position: 1
        
  - id: "#reference"
    type: File
    inputBinding:
      position: 2
      prefix: "REFERENCE="
      secondaryFiles:
        - ".amb"
        - ".ann"
        - ".bwt"
        - ".pac"
        - ".sa"
        
  - id: "#outputFileName"
    type: string
    inputBinding:
      position: 3
      prefix: "OUTPUT="

  - id: "#genomeAssembly"
    type: string
    default: "null"
    inputBinding:
      position: 4
      prefix: "GENOME_ASSEMBLY="

  - id: "#uri"
    type: string
    default: "null"
    inputBinding:
      position: 5
      prefix: "URI="

  - id: "#species"
    type: string
    default: "null"
    inputBinding:
      position: 6
      prefix: "SPECIES="  

  - id: "#truncateName"
    type: ["null", string]
    default: "true"
    inputBinding:
      position: 7
      prefix: "TRUNCATE_NAMES_AT_WHITESPACE="
      
  - id: "#numSequences"
    type: int
    default: 2147483647
    inputBinding:
      position: 8
      prefix: "NUM_SEQUENCES="
      
outputs:
  - id: "#createDict_output"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputFileName
      
baseCommand: []
      
      


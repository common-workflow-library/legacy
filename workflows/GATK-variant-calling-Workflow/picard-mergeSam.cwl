#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool

description: |
  Merges multiple SAM/BAM files into one file
  Options: INPUT File SAM or BAM input file Default value null. This option must be specified at least 1 times
    OUTPUT File SAM or BAM file to write merged result to Required
    SORT_ORDER SortOrder Sort order of output file Default value coordinate. This option can be set to null to clear the default value. Possible values {unsorted, queryname, coordinate, duplicate}
    ASSUME_SORTED Boolean If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    MERGE_SEQUENCE_DICTIONARIES Boolean Merge the sequence dictionaries Default value false. This option can be set to null to clear the default value. Possible values {true, false}
    USE_THREADING Boolean Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    COMMENT String Comment(s) to include in the merged output files header. Default value null. This option may be specified 0 or more times
   
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: picard-docker.cwl
  
inputs:

  - id: "#mergeSam"
    type: string
    default: "MergeSamFiles"
    inputBinding:
      position: 1

  - id: "#outputFileName_mergedSam"
    type: string
    inputBinding:
      prefix: "OUTPUT="

  - id: "#inputFileName_mergedSam"
    type:
      type: array
      items: File
      inputBinding: { prefix: "INPUT=" } 

  - id: "#sortorder"
    type: string
    default: "coordinate"
    inputBinding:
      prefix: "SORT_ORDER="

  - id: "#readSorted"
    type: ["null", string]
    inputBinding:
      prefix: "ASSUME_SORTED="

  - id: "#mergeSequenceDictionaries"
    type: ["null", string]
    inputBinding:
      prefix: "MERGE_SEQUENCE_DICTIONARIES="
      
  - id: "#useThreading"
    type: ["null", string]
    inputBinding:
      prefix: "USE_THREADING="
      
  - id: "#comment"
    type: string
    inputBinding:
      prefix: "COMMENT="

  - id: "#tmpdir"
    type: string
    inputBinding:
      prefix: "TMP_DIR="

outputs:
  - id: "#mergeSam_output"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputFileName_mergedSam
      
baseCommand: []

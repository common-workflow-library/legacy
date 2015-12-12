#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool

description: |
  Sorts the input SAM or BAM. Input and output formats are determined by file extension.
  Options: INPUT File The BAM or SAM file to sort. Required
    OUTPUT File The sorted BAM or SAM output file. Required
    SORT_ORDER SortOrder Sort order of output file Required. Possible values {unsorted, queryname, coordinate, duplicate}
   
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: picard-docker.cwl
  
inputs:

  - id: "#sortSam"
    type: string
    default: "SortSam"
    inputBinding:
      position: 1
      
  - id: "#inputFileName_sortSam"
    type: File
    inputBinding:
      position: 2
      prefix: "INPUT="
      
  - id: "#outputFileName_sortSam"
    type: string
    inputBinding:
      position: 3
      prefix: "OUTPUT="

  - id: "#SO-coordinate"
    type: string
    default: "coordinate"
    inputBinding:
      position: 4
      prefix: "SORT_ORDER="

  - id: "#tmpdir"
    type: string
    inputBinding:
      position: 5
      prefix: "TMP_DIR="

  - id: "#createIndex"
    type: ["null", string]
    default: 'true'
    inputBinding:
      position: 6
      prefix: "CREATE_INDEX="

      
outputs:
  - id: "#sortSam_output"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputFileName_sortSam
      
baseCommand: []

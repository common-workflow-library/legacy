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
  
inputs:
  - id: "#java_arg"
    type: string
    default: "-Xmx2g"
    inputBinding:
      position: 1
      
  - id: "#jar_file"
    type: File
    inputBinding:
      position: 2
      prefix: "-jar"

  - id: "#sortSam"
    type: string
    default: "SortSam"
    inputBinding:
      position: 3
      
  - id: "#inputFileName_mergedSam"
    type: File
    inputBinding:
      position: 4
      prefix: "INPUT="
      
  - id: "#outputFileName_sortSam"
    type: string
    inputBinding:
      position: 5
      prefix: "OUTPUT="

  - id: "#SO-coordinate"
    type: string
    default: "coordinate"
    inputBinding:
      position: 6
      prefix: "SO="
      
outputs:
  - id: "#sortSam_output"
    type: File
    outputBinding: {glob: "sortedSam.bam"}
      
baseCommand: ["java"]

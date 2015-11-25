#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool

description: |
  Generates a BAM index (.bai) file.
  Options: 
    INPUT String A BAM file or URL to process. Must be sorted in coordinate order. Required
    OUTPUT File The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai. If INPUT is a URL and OUTPUT is unspecified, defaults to a file in the current directory. Default value null
   
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

  - id: "#BuildBamIndex"
    type: string
    default: "BuildBamIndex"
    inputBinding:
      position: 3
      
  - id: "#outputFileName_BuildBamIndex"
    type: string
    inputBinding:
      position: 4
      prefix: "OUTPUT="
      
  - id: "#inputFileName_BuildBamIndex"
    type: File
    inputBinding:
      position: 5
      prefix: "INPUT="
      
outputs:
  - id: "#markDups_output"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputFileName_BuildBamIndex
      
baseCommand: ["java"]

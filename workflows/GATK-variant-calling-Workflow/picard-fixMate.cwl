#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool

description: |
  Ensure that all mate-pair information is in sync between each read and its mate pair. If no OUTPUT file is supplied then the output is written to a temporary file and then copied over the INPUT file. Reads marked with the secondary alignment flag are written to the output file unchanged.  
  Options:
    INPUT File The input file to fix. Default value null. This option may be specified 0 or more times.
    OUTPUT File The output file to write to. If no output file is supplied, the input file is overwritten. Default value null.
    SORT_ORDER SortOrder Optional sort order if the OUTPUT file should be sorted differently than the INPUT file. Default value null. Possible values: {unsorted, queryname, coordinate, duplicate}
    ASSUME_SORTED Boolean If true, assume that the input file is queryname sorted, even if the header says otherwise. Default value: false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    ADD_MATE_CIGAR Boolean Adds the mate CIGAR tag (MC) if true, does not if false. Default value true. This option can be set to 'null' to clear the default value. Possible values {true, false}
    IGNORE_MISSING_MATES  Boolean If true, ignore missing mates, otherwise will throw an exception when missing mates are found. Default value: true. This option can be set to 'null' to clear the default value. Possible values {true, false}
    
    
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: picard-docker.cwl

inputs:

  - id: "#FixMateInformation"
    type: string
    default: "FixMateInformation"
    inputBinding:
      position: 1

  - id: "#inputFileName_fixMate"
    type: File
    inputBinding:
      position: 2
      prefix: "INPUT="
        
  - id: "#outputFileName_fixMate"
    type: string
    inputBinding:
      position: 3
      prefix: "OUTPUT="

  - id: "#sortOrder"
    type: string
    default: "coordinate"
    inputBinding:
      position: 4
      prefix: "SORT_ORDER="
      
  - id: "#createIndex"
    type: ["null", string]
    default: 'true'
    inputBinding:
      position: 5
      prefix: "CREATE_INDEX="

  - id: "#tmpdir"
    type: string
    inputBinding:
      position: 6
      prefix: "TMP_DIR="
      
  - id: "#readSorted"
    type: ["null", string]
    inputBinding:
      position: 7
      prefix: "ASSUME_SORTED="
      
  - id: "#mateCigar"
    type: ["null", string]
    inputBinding:
      position: 8
      prefix: "ADD_MATE_CIGAR="
      
  - id: "#ignoreMissingMates"
    type: ["null", string]
    inputBinding:
      position: 9
      prefix: "IGNORE_MISSING_MATES="
      
outputs:
  - id: "#fixMate_output"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputFileName_fixMate
    
baseCommand: []





 

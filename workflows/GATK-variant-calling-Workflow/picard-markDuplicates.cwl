#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool

description: |
  Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged
  Options:
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP INT Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value 8000. This option can be set to 'null' to clear the default value
    SORTING_COLLECTION_SIZE_RATIO DOUBLE This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections. If you are running out of memory, try reducing this number. Default value 0.25. This option can be set to 'null' to clear the default value
    BARCODE_TAG STR Barcode SAM tag (ex. BC for 10X Genomics) Default value null
    READ_ONE_BARCODE_TAG STR Read one barcode SAM tag (ex. BX for 10X Genomics) Default value null
    READ_TWO_BARCODE_TAG STR Read two barcode SAM tag (ex. BX for 10X Genomics) Default value null
    INPUT STR One or more input SAM or BAM files to analyze. Must be coordinate sorted. Default value null. This option may be specified 0 or more times
    OUTPUT FILE The output file to write marked records to Required
    METRICS_FILE FILE File to write duplication metrics to Required
    PROGRAM_RECORD_ID STR The program record ID for the @PG record(s) created by this program. Set to null to disable PG record creation. This string may have a suffix appended to avoid collision with other program record IDs. Default value MarkDuplicates. This option can be set to 'null' to clear the default value
    PROGRAM_GROUP_VERSION  STR Value of VN tag of PG record to be created. If not specified, the version will be detected automatically. Default value null
    PROGRAM_GROUP_COMMAND_LINE STR Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically. Default value null
    PROGRAM_GROUP_NAME STR Value of PN tag of PG record to be created. Default value MarkDuplicates. This option can be set to 'null' to clear the default value
    COMMENT STR Comment(s) to include in the output files header. Default value null. This option may be specified 0 or more times
    REMOVE_DUPLICATES BOOLEAN If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    ASSUME_SORTED BOOLEAN If true, assume that the input file is coordinate sorted even if the header says otherwise. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    DUPLICATE_SCORING_STRATEGY SCORINGSTRATEGY The scoring strategy for choosing the non-duplicate among candidates. Default value SUM_OF_BASE_QUALITIES. This option can be set to 'null' to clear the default value. Possible values {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH}
    READ_NAME_REGEX STR Regular expression that can be used to parse read names in the incoming SAM file. Read names are parsed to extract three variables tile/region, x coordinate and y coordinate. These values are used to estimate the rate of optical duplication in order to give a more accurate estimated library size. Set this option to null to disable optical duplicate detection. The regular expression should contain three capture groups for the three variables, in order. It must match the entire read name. Note that if the default regex is specified, a regex match is not actually done, but instead the read name is split on colon character. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values. Default value [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*. This option can be set to 'null' to clear the default value
    OPTICAL_DUPLICATE_PIXEL_DISTANCE INT The maximum offset between two duplicte clusters in order to consider them optical duplicates. This should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal. Default value 100. This option can be set to 'null' to clear the default value

requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: picard-docker.cwl

inputs:

  - id: "#MarkDuplicates"
    type: string
    default: "MarkDuplicates"
    inputBinding:
      position: 1

  - id: "#inputFileName_markDups"
    type: File
    inputBinding:
      position: 2
      prefix: "INPUT="
      secondaryFiles:
        - "^.bai"
 
  - id: "#outputFileName_markDups"
    type: string
    inputBinding:
      position: 3
      prefix: "OUTPUT="

  - id: "#metricsFile"
    type: string
    inputBinding:
      position: 4
      prefix: "METRICS_FILE="

  - id: "#readSorted"
    type: ["null", string]
    inputBinding:
      position: 5
      prefix: "ASSUME_SORTED="

  - id: "#removeDuplicates"
    type: ["null", string]
    inputBinding:
      position: 6
      prefix: "REMOVE_DUPLICATES="
      
  - id: "#maxFileHandles"
    type: int
    default: 8000
    inputBinding:
      position: 7
      prefix: "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP="
      
  - id: "#sortRatio"
    type: double
    default: 0.25
    inputBinding:
      position: 8
      prefix: "SORTING_COLLECTION_SIZE_RATIO="
      
  - id: "#barcodeTag"
    type: ["null", string]
    inputBinding:
      position: 9
      prefix: "BARCODE_TAG="
      
  - id: "#readOneBarcodeTag"
    type: ["null", string]
    inputBinding:
      position: 10
      prefix: "READ_ONE_BARCODE_TAG="
      
  - id: "#readTwoBarcodeTag"
    type: ["null", string]
    inputBinding:
      position: 11
      prefix: "READ_TWO_BARCODE_TAG="
      
  - id: "#recordId"
    type: string
    default: 'MarkDuplicates'
    inputBinding:
      position: 12
      prefix: "PROGRAM_RECORD_ID="
      
  - id: "#groupVersion"
    type: ["null", string]
    inputBinding:
      position: 13
      prefix: "PROGRAM_GROUP_VERSION="
      
  - id: "#groupCommandLine"
    type: ["null", string]
    inputBinding:
      position: 14
      prefix: "PROGRAM_GROUP_COMMAND_LINE="
      
  - id: "#groupCommandName"
    type: string
    default: 'MarkDuplicates'
    inputBinding:
      position: 15
      prefix: "PROGRAM_GROUP_NAME="
      
  - id: "#comment"
    type: ["null", string]
    inputBinding:
      position: 16
      prefix: "COMMENT="
      
  - id: "#regularExpression"
    type: string
    default: '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*'
    inputBinding:
      position: 17
      prefix: "READ_NAME_REGEX="
      
  - id: "#pixelDistance"
    type: int
    default: 100
    inputBinding:
      position: 18
      prefix: "OPTICAL_DUPLICATE_PIXEL_DISTANCE="
      
  - id: "#createIndex"
    type: ["null", string]
    default: 'true'
    inputBinding:
      position: 19
      prefix: "CREATE_INDEX="

  - id: "#tmpdir"
    type: string
    inputBinding:
      position: 20
      prefix: "TMP_DIR="
      
outputs:
  - id: "#markDups_output"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputFileName_markDups
    
baseCommand: []





 

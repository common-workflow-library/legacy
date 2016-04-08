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
    http://broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex
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
  picard-BuildBamIndex.cwl is developed for CWL consortium
    Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged

doap:name: "picard-MarkDuplicates.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/picard-MarkDuplicates.cwl"

dct:creator:
- class: foaf:Organization
  foaf:name: "UNIVERSITY OF MELBOURNE"
  foaf:member:
  - class: foaf:Person
    id: "farahk@student.unimelb.edu.au"
    foaf:mbox: "mailto:farahk@student.unimelb.edu.au"
  - class: foaf:Person
    id: "skanwal@student.unimelb.edu.au"
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
- $import: picard-docker.yml
- class: InlineJavascriptRequirement

inputs:

  - id: "java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding:
      position: 1
      
  - id: "#inputFileName_markDups"
    description: One or more input SAM or BAM files to analyze. Must be coordinate sorted. Default value null. This option may be specified 0 or more times
    type:
      - "null"
      - type: array
        items: File
        inputBinding: 
          prefix: "INPUT="
        secondaryFiles:
          - "^.bai"
    inputBinding:
      position: 4
 
  - id: "#outputFileName_markDups"
    type: string
    description: The output file to write marked records to Required
    inputBinding:
      position: 5
      prefix: "OUTPUT="

  - id: "#metricsFile"
    type: ["null", string]
    description: File to write duplication metrics to Required
    inputBinding:
      position: 6
      prefix: "METRICS_FILE="


  - id: "#removeDuplicates"
    type: ["null", boolean]
    description: If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    inputBinding:
      position: 7
      prefix: "REMOVE_DUPLICATES="
      
  - id: "#maxFileHandles"
    type: ["null", int]
    description: Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value 8000. This option can be set to 'null' to clear the default value
    inputBinding:
      position: 8
      prefix: "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP="
      
  - id: "#sortRatio"
    type: ["null", double]
    description: This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections. If you are running out of memory, try reducing this number. Default value 0.25. This option can be set to 'null' to clear the default value
    inputBinding:
      position: 9
      prefix: "SORTING_COLLECTION_SIZE_RATIO="
      
  - id: "#barcodeTag"
    type: ["null", string]
    description: Barcode SAM tag (ex. BC for 10X Genomics) Default value null
    inputBinding:
      position: 10
      prefix: "BARCODE_TAG="
      
  - id: "#readOneBarcodeTag"
    type: ["null", string]
    description: Read one barcode SAM tag (ex. BX for 10X Genomics) Default value null
    inputBinding:
      position: 11
      prefix: "READ_ONE_BARCODE_TAG="
      
  - id: "#readTwoBarcodeTag"
    type: ["null", string]
    description: Read two barcode SAM tag (ex. BX for 10X Genomics) Default value null
    inputBinding:
      position: 12
      prefix: "READ_TWO_BARCODE_TAG="
      
  - id: "#recordId"
    type: ["null", string]
    description: The program record ID for the @PG record(s) created by this program. Set to null to disable PG record creation. This string may have a suffix appended to avoid collision with other program record IDs. Default value MarkDuplicates. This option can be set to 'null' to clear the default value
    inputBinding:
      position: 13
      prefix: "PROGRAM_RECORD_ID="
      
  - id: "#groupVersion"
    type: ["null", string]
    description: Value of VN tag of PG record to be created. If not specified, the version will be detected automatically. Default value null
    inputBinding:
      position: 14
      prefix: "PROGRAM_GROUP_VERSION="
      
  - id: "#groupCommandLine"
    type: ["null", string]
    description: Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically. Default value null
    inputBinding:
      position: 15
      prefix: "PROGRAM_GROUP_COMMAND_LINE="
      
  - id: "#groupCommandName"
    type: ["null", string]
    description: Value of PN tag of PG record to be created. Default value MarkDuplicates. This option can be set to 'null' to clear the default value
    inputBinding:
      position: 16
      prefix: "PROGRAM_GROUP_NAME="
      
  - id: "#comment"
    description: Comment(s) to include in the output files header. Default value null. This option may be specified 0 or more times
    type:
      - "null"
      - type: array
        items: File
        inputBinding: 
          prefix: "COMMENT="
    inputBinding:
      position: 17
      
  - id: "#regularExpression"
    type: ["null", string]
    description: Regular expression that can be used to parse read names in the incoming SAM file. Read names are parsed to extract three variables tile/region, x coordinate and y coordinate. These values are used to estimate the rate of optical duplication in order to give a more accurate estimated library size. Set this option to null to disable optical duplicate detection. The regular expression should contain three capture groups for the three variables, in order. It must match the entire read name. Note that if the default regex is specified, a regex match is not actually done, but instead the read name is split on colon character. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values. Default value [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*. This option can be set to 'null' to clear the default value
    inputBinding:
      position: 18
      prefix: "READ_NAME_REGEX="
      
  - id: "#pixelDistance"
    type: ["null", int]
    description: The maximum offset between two duplicte clusters in order to consider them optical duplicates. This should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal. Default value 100. This option can be set to 'null' to clear the default value
    inputBinding:
      position: 19
      prefix: "OPTICAL_DUPLICATE_PIXEL_DISTANCE="
      
  - id: "#createIndex"
    type: ["null", string]
    default: "true"
    description: Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    inputBinding:
      position: 20
      prefix: "CREATE_INDEX="

  - id: "#tmpdir"
    type: ["null", string]
    description: Default value null. This option may be specified 0 or more times.
    inputBinding:
      position: 21
      prefix: "TMP_DIR="

  - id: "#readSorted"
    type: ["null", boolean]
    description: If true, assume that the input file is coordinate sorted even if the header says otherwise. Default value false. This option can be set to 'null' to clear the default value. Possible values {true, false}
    inputBinding:
      position: 22
      prefix: "ASSUME_SORTED="
      
outputs:
  - id: "markDups_output"
    type: File
    outputBinding: 
      glob: $(inputs.outputFileName_markDups)
    
baseCommand: ["java"]

arguments: 

  - valueFrom: "/usr/local/bin/picard.jar"
    position: 2
    prefix: "-jar"
  - valueFrom: "MarkDuplicates"
    position: 3

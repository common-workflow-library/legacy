#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

#!/usr/bin/env cwl-runner
class: CommandLineTool

description: |
  Prints all reads that have a mapping quality above zero
  Usage: java -Xmx4g -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta -I input1.bam -I input2.bam -o output.bam --read_filter MappingQualityZero
  Options: 
   -T PrintReads Tool name to be executed
   -R reference.fasta Reference sequence in fasta format
   -I input.bam bam file produce from indelRealigner
   -BQSR recalibration_report.grp the recalibration table data in recalibration_report.grp produced by BaseRecalibration to recalibrate the quality scores in input.bam
   -o output.bam Output file name

requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: gatk-docker.cwl

inputs:

  - id: "#java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding: 
      position: 1

  - id: "#jar_file"
    type: string
    inputBinding: { position: 2, prefix: "-jar" }
    description: GATK jar file

     
  - id: "#printReads"
    type: string
    default: "PrintReads"
    inputBinding: { position: 3, prefix: "-T" }
    description: tool used for this step from GATK jar
     
  - id: "#reference"
    type: File
    inputBinding:
      position: 4
      prefix: "-R"
      secondaryFiles:
        - ".amb"
        - ".ann"
        - ".bwt"
        - ".pac"
        - ".rbwt"
        - ".rpac"
        - ".rsa"
        - ".sa"
        - ".fai"
        - "^.dict"

  - id: "#inputBam_printReads"
    type: File
    description: bam file produced after indelRealigner
    inputBinding:
      position: 5 
      prefix: "-I"
      secondaryFiles:
        - "^.bai"
      
  - id: "#input_baseRecalibrator"
    type: File
    description: the recalibration table produced by BaseRecalibration
    inputBinding:
      position: 5 
      prefix: "-BQSR"
      
  - id: "#outputfile_printReads"
    type: string
    description: name of the output file from indelRealigner
    inputBinding:
      position: 9
      prefix: "-o"
    
outputs:
  - id: "#output_printReads"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputfile_printReads

baseCommand: ["java"]


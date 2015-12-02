#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

#!/usr/bin/env cwl-runner
class: CommandLineTool

description: |
  Generate base recalibration table to compensate for systematic errors in basecalling confidences
  Options:
   -T BaseRecalibrator Tool name to be executed
   -R reference.fasta Reference sequence in fasta format
   -I my_reads.bam bam file produced from indelRealigner
   -knownSites latest_dbsnp.vcf set of known indels
   -o recal_data.table Output file name
   
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

     
  - id: "#BaseRecalibrator"
    type: string
    default: "BaseRecalibrator"
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

  - id: "#inputBam_BaseRecalibrator"
    type: File
    description: bam file produced after indelRealigner
    inputBinding:
      position: 5 
      prefix: "-I"
      secondaryFiles:
        - "^.bai"      
     
  - id: "#known"  
    type:
      type: array
      items: File
      inputBinding: { prefix: "-knownSites" }
    inputBinding: 
      position: 6
      
  - id: "#outputfile_BaseRecalibrator"
    type: string
    description: name of the output file from baseRecalibrator
    inputBinding:
      position: 7
      prefix: "-o"
    
outputs:
  - id: "#output_baseRecalibrator"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputfile_BaseRecalibrator

baseCommand: ["java"]


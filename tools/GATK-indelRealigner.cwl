#!/usr/bin/env cwl-runner
class: CommandLineTool

description: |
  Perform local realignment of reads around indels
  Options: 
   -T IndelRealigner Tool name to be executed
   -R reference.fasta Reference sequence in fasta format
   -I input.bam bam file produced from realignerTargetCreator
   -known indels.vcf set of known indels 
   -targetIntervals intervalListFromRTC.intervals Restrits alignment to provided intervals. Output dataset from realignerTargetCreator
   -o realignedBam.bam Output file name

requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: gatk-docker.cwl
arguments:
  - valueFrom: "/home/biodocker/bin/gatk/target/GenomeAnalysisTK.jar"
    position: 2
    prefix: "-jar"
    
inputs:

  - id: "#java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding: 
      position: 1
     
  - id: "#IndelRealigner"
    type: string
    default: "IndelRealigner"
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
        
  - id: "#inputBam_realign"
    type: File
    description: bam file produced after markDups execution
    inputBinding:
      position: 5 
      prefix: "-I"
      secondaryFiles:
        - ".bai"
        
  - id: "#known"  
    type:
      type: array
      items: File
      inputBinding: { prefix: "-known" }
    inputBinding: 
      position: 6
      
  - id: "#intervals"
    type: File
    description: list of intervals created by realignerTargetCreataor
    inputBinding:
      position: 7
      prefix: "-targetIntervals"
      
  - id: "#outputfile_indelRealigner"
    type: string
    description: name of the output file from indelRealigner
    inputBinding:
      position: 9
      prefix: "-o"
    
outputs:
  - id: "#output_indelRealigner"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputfile_indelRealigner

baseCommand: ["java"]

#!/usr/bin/env cwl-runner

class: CommandLineTool
requirements:
  - import: node-engine.cwl
description: |
"This is a tool wrapper for GATK tool called realignerTargetCreator. it accepts 4 inputs and produces a file containing list of target intervals to pass to the IndelRealigner."
	Usage: java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference fasta -I Input.bam --known indels.vcf -o forIndelRealigner.intervals"

inputs:
  - id: "#toolname"
    type: File
    inputBinding: { position: 1, prefix: "-T" }   
  - id: "#reference"
    type: File
    inputBinding: { position: 2, prefix: "-R" } #where to specify the file type? 

  - id: "#inputfile"
    type: File
    inputBinding: { position: 3, prefix: "-I" }
    
  - id: "#known"
    type: File
    inputBinding: { position: 4, prefix: "--known" }
   
outputs:
  - id: "#forIndelRealigner"
    type: "File"
    outputBinding: { "glob": "forIndelRealigner.intervals", prefix: "-o" }

baseCommand: ["java"]

stdout: "forIndelRealigner.intervals"

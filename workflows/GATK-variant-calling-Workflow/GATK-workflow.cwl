#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/

class: Workflow

inputs:
  - id: "#reference"
    type: File
    description: reference human genome file 
    
  - id: "#FwdReadFile"
    type: File
    description: file containing the first end of paired end reads in fastq format

  - id: "#RevReadFile"
    type: File
    description: file containing the second end of paired end reads in fastq format

  - id: "#output_name"
    type: string
    description: name of bwa-mem output file

  - id: "#outputFileName"
    type: string
    description: output file name for picard create dictionary command from picard toolkit
    
  - id: "#outputFileName_mergedSam"
    type: string
    description: output file name generated as a result of mergesam command from picard toolkit

  - id: "#outputFileName_sortSam"
    type: string
    description: output file name generated as a result of sortSam command from picard toolkit

  - id: "#outputFileName_fixMate"
    type: string
    description: output file name generated as a result of fixMate command from picard toolkit

  - id: "#outputFileName_markDups"
    type: string
    description: output file name generated as a result of Markduplicates command from picard toolkit
    
  - id: "#metricsFile"
    type: string
    description: metric file generated by MarkDupicates command listing duplicates

  - id: "#readSorted"
    type: string
    description: set to be true showing that reads are sorted already

  - id: "#removeDuplicates"
    type: string
    description: set to be true
 
  - id: "#outputfile_realignTarget"
    type: string
    description: name of realignTargetCreator output file
    
  - id: "#createIndex"
    type: string
    description: set to be true to create .bai file
    
  - id: "#tmpdir"
    type: string
    description: temporary directory for markdups


  - id: "#known"
    type:
      - "null"
      - type: array
        items: File
    description: array of known variant files for realign target creator

  - id: "#outputfile_indelRealigner"
    type: string
    description: name of indelRealigner output file

  - id: "#outputfile_BaseRecalibrator"
    type: string
    description: name of BaseRecalibrator output file

  - id: "#outputfile_printReads"
    type: string
    description: name of PrintReads command output file
    
  - id: "#outputfile_HaplotypeCaller"
    type: string
    description: name of Haplotype caller command output file

  - id: "#dbsnp"
    type: File
    description: vcf file containing SNP variations used for Haplotype caller


outputs:
  - id: "#bwamem_output"
    type: File
    source: "#bwa-mem.readfileSam_output"

  - id: "#createDict_output"
    type: File
    source: "#create-dict.createDict_output"

  - id: "#mergedSam_output"
    type: File
    source: "#mergeSam.mergeSam_output"
    
  - id: "#sortSam_output"
    type: File
    source: "#sortSam.sortSam_output"

  - id: "#fixMate_output"
    type: File
    source: "#fixMate.fixMate_output"

  - id: "#markDups_output"
    type: File
    source: "#markDups.markDups_output"
    
  - id: "#output_realignTarget"
    type: File
    source: "#RealignTarget.output_realignTarget" 

  - id: "#output_indelRealigner"
    type: File
    source: "#IndelRealigner.output_indelRealigner"
    
  - id: "#outputfile_baseRecalibrator"
    type: File
    source: "#BaseRecalibrator.output_baseRecalibrator"

  - id: "#output_printReads"
    type: File
    source: "#PrintReads.output_printReads"
    
  - id: "#output_HaplotypeCaller"
    type: File
    source: "#HaplotypeCaller.output_HaplotypeCaller"



requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl

steps:
  - id: "#bwa-mem"
    run: { import: bwa-mem.cwl }
    inputs:

      - { id: "#bwa-mem.reference", source: "#reference" }
      - { id: "#bwa-mem.FwdReadFile", source: "#FwdReadFile" }
      - { id: "#bwa-mem.RevReadFile", source: "#RevReadFile" }
      - { id: "#bwa-mem.output_name", source: "#output_name" }
    outputs:
      - { id: "#bwa-mem.readfileSam_output" }
  
  - id: "#create-dict"
    run: { import: createDict.cwl }
    inputs:

      - { id: "#create-dict.reference", source: "#reference" }
      - { id: "#create-dict.outputFileName", source: "#outputFileName" }
      
    outputs:
      - { id: "#create-dict.createDict_output" }

  - id: "#mergeSam"
    run: { import: mergeSam.cwl }
    inputs:
      - { id: "#mergeSam.outputFileName_mergedSam", source: "#outputFileName_mergedSam" }
      - { id: "#mergeSam.inputFileName_mergedSam", source: ["#bwa-mem.readfileSam_output", "#create-dict.createDict_output" ]}      
    outputs:
      - { id: "#mergeSam.mergeSam_output" }

  - id: "#sortSam"
    run: { import: sortSam.cwl }
    inputs:
      - { id: "#sortSam.outputFileName_sortSam", source: "#outputFileName_sortSam" }
      - { id: "#sortSam.inputFileName_sortSam", source: "#mergeSam.mergeSam_output" }
      - { id: "#sortSam.createIndex", source: "#createIndex" } 
      - { id: "#sortSam.tmpdir", source: "#tmpdir" } 
      
    outputs:
      - { id: "#sortSam.sortSam_output" }

  - id: "#fixMate"
    run: { import: fixMate.cwl }
    inputs:
      - { id: "#fixMate.outputFileName_fixMate", source: "#outputFileName_fixMate" }
      - { id: "#fixMate.inputFileName_fixMate", source: "#sortSam.sortSam_output" }
      - { id: "#fixMate.tmpdir", source: "#tmpdir" } 
      
    outputs:
      - { id: "#fixMate.fixMate_output" }
      
  - id: "#markDups"
    run: { import: markDups.cwl }
    inputs:
      - { id: "#markDups.outputFileName_markDups", source: "#outputFileName_markDups" }
      - { id: "#markDups.inputFileName_markDups", source: "#fixMate.fixMate_output" }
      - { id: "#markDups.metricsFile", source: "#metricsFile" }
      - { id: "#markDups.readSorted", source: "#readSorted" }
      - { id: "#markDups.removeDuplicates", source: "#removeDuplicates" } 
      - { id: "#markDups.createIndex", source: "#createIndex" } 
      - { id: "#markDups.tmpdir", source: "#tmpdir" } 
      
    outputs:
      - { id: "#markDups.markDups_output" }

  - id: "#RealignTarget"
    run: { import: realignTargetCreator.cwl }
    inputs:
      - { id: "#RealignTarget.outputfile_realignTarget", source: "#outputfile_realignTarget" }
      - { id: "#RealignTarget.inputBam_realign", source: "#markDups.markDups_output" }
      - { id: "#RealignTarget.reference", source: "#reference" }
      - { id: "#RealignTarget.known", source: "#known" }
      
    outputs:
      - { id: "#RealignTarget.output_realignTarget" }

  - id: "#IndelRealigner"
    run: { import: indelRealigner.cwl }
    inputs:
      - { id: "#IndelRealigner.outputfile_indelRealigner", source: "#outputfile_indelRealigner" }
      - { id: "#IndelRealigner.inputBam_realign", source: "#markDups.markDups_output" }
      - { id: "#IndelRealigner.intervals", source: "#RealignTarget.output_realignTarget" }
      - { id: "#IndelRealigner.reference", source: "#reference" }
      - { id: "#IndelRealigner.known", source: "#known" }
      
    outputs:
      - { id: "#IndelRealigner.output_indelRealigner" }

  - id: "#BaseRecalibrator"
    run: { import: BaseRecalibrator.cwl }
    inputs:
      - { id: "#BaseRecalibrator.outputfile_BaseRecalibrator", source: "#outputfile_BaseRecalibrator" }
      - { id: "#BaseRecalibrator.inputBam_BaseRecalibrator", source: "#IndelRealigner.output_indelRealigner" }
      - { id: "#BaseRecalibrator.reference", source: "#reference" }
      - { id: "#BaseRecalibrator.known", source: "#known" }
      
    outputs:
      - { id: "#BaseRecalibrator.output_baseRecalibrator" }

  - id: "#PrintReads"
    run: { import: printReads.cwl }
    inputs:
      - { id: "#PrintReads.outputfile_printReads", source: "#outputfile_printReads" }
      - { id: "#PrintReads.inputBam_printReads", source: "#IndelRealigner.output_indelRealigner" }
      - { id: "#PrintReads.reference", source: "#reference" }
      - { id: "#PrintReads.input_baseRecalibrator", source: "#BaseRecalibrator.output_baseRecalibrator" }
    outputs:
      - { id: "#PrintReads.output_printReads" }
      
  - id: "#HaplotypeCaller"
    run: { import: haplotypeCaller.cwl }
    inputs:
      - { id: "#HaplotypeCaller.outputfile_HaplotypeCaller", source: "#outputfile_HaplotypeCaller" }
      - { id: "#HaplotypeCaller.inputBam_HaplotypeCaller", source: "#PrintReads.output_printReads" }
      - { id: "#HaplotypeCaller.reference", source: "#reference" }
      - { id: "#HaplotypeCaller.dbsnp", source: "#dbsnp" }
    outputs:
      - { id: "#HaplotypeCaller.output_HaplotypeCaller" }


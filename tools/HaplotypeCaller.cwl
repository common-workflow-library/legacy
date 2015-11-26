#!/usr/bin/env cwl-runner
#
# Authors: farahk@student.unimelb.edu.au and skanwal@student.unimelb.edu.au UNIVERSITY OF MELBOURNE
# Developed for CWL consortium http://commonwl.org/
class: CommandLineTool

description: |
  Call germline SNPs and indels via local re-assembly of haplotypes
  Options:
    -R reference.fasta Reference sequence in fasta format
    -T HaplotypeCaller Tool name to be executed
    -I sample.bam bam file from printReads
    [--dbsnp dbSNP.vcf] latest_dbsnp.vcf set of known indels
    [-stand_call_conf 30] The minimum phred-scaled confidence threshold at which variants should be called - optional
    [-stand_emit_conf 10] The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold) - optonal
    [-L targets.interval_list] target bed file - if any
    -o output.raw.snps.indels.vcf  Output file from haplotypr caller
   
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl

inputs:

  - id: "#java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding: 
      position: 1

  - id: "#jar_file"
    type: File
    inputBinding: { position: 2, prefix: "-jar" }
    description: GATK jar file

     
  - id: "#HaplotypeCaller"
    type: string
    default: "HaplotypeCaller"
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

  - id: "#inputBam_HaplotypeCaller"
    type: File
    description: bam file produced after printReads
    inputBinding:
      position: 5 
      prefix: "-I"
    
  - id: "#dbsnp"  
    type: File
    description: latest_dbsnp.vcf set of known indels
    inputBinding: 
      position: 6
      prefix: "--dbsnp"
      
  - id: "#outputfile_HaplotypeCaller"
    type: string
    description: name of the output file from HaplotypeCaller
    inputBinding:
      position: 7
      prefix: "-o"
    
outputs:
  - id: "#output_HaplotypeCaller"
    type: File
    outputBinding: 
      glob:
        engine: cwl:JsonPointer
        script: /job/outputfile_HaplotypeCaller

baseCommand: ["java"]


#!/usr/bin/env cwl-runner

class: CommandLineTool

description: |
  Usage:  samtools rmdup [-sS] <input.srt.bam> <output.bam>

  Option: -s    rmdup for SE reads
          -S    treat PE reads as SE in rmdup (force -s)


requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: samtools-docker.cwl

inputs:
  - id: "#input"
    type: File
    description: |
      Input bam file.
    inputBinding:
      position: 2

  - id: "#output_name"
    type: string
    inputBinding:
      position: 3

  - id: "#single_end"
    type: boolean
    default: false
    description: |
      rmdup for SE reads

  - id: "#pairend_as_se"
    type: boolean
    default: false
    description: |
      treat PE reads as SE in rmdup (force -s)

outputs:
  - id: "#rmdup"
    type: File
    description: "File with removed duplicates"
    outputBinding:
      glob:
        engine: cwl:JsonPointer
        script: /job/output_name

baseCommand: ["samtools", "rmdup"]

#!/usr/bin/env cwl-runner

class: CommandLineTool

description: |
  Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
  Options:
    -b       Generate BAI-format index for BAM files [default]
    -c       Generate CSI-format index for BAM files
    -m INT   Set minimum interval size for CSI indices to 2^INT [14]

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
      position: 4

  - id: "#fakeoutput"
    type: string
    default: ""
    inputBinding:
      position: 6
      valueFrom:
        engine: node-engine.cwl
        script: |
          {
            var ext=$job['bai']?'.bai':$job['csi']?'.csi':'.bai';
            return $job['input'].path.split('/').slice(-1)[0]+ext;
          }

  - id: "#bai"
    type: boolean
    default: false
    description: |
      Generate BAI-format index for BAM files [default]

  - id: "#csi"
    type: boolean
    default: false
    description: |
      Generate CSI-format index for BAM files

  - id: "#interval"
    type: ["null", int]
    description: |
      Generate CSI-format index for BAM files
    inputBinding:
      position: 1
      prefix: "-m"

outputs:
  - id: "#sorted"
    type: File
    description: "The sorted file"
    outputBinding:
      glob:
        engine: node-engine.cwl
        script: |
          {
            var ext=$job['bai']?'.bai':$job['csi']?'.csi':'.bai';
            return $job['input'].path.split('/').slice(-1)[0]+ext;
          }

baseCommand: ["samtools", "index"]

arguments:
  - valueFrom:
      engine: node-engine.cwl
      script: |
        $job['bai']?'-b':$job['csi']?'-c':[]

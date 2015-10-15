#!/usr/bin/env cwl-runner

class: CommandLineTool

description: |
  Usage: bwa samse [-n max_occ] [-f out.sam] [-r RG_line] <prefix> <in.sai> <in.fq>
requirements:
  - import: node-engine.cwl
  - import: envvar-global.cwl
  - import: bwa-docker.cwl

inputs:
  - id: "#prefix"
    type: File
    inputBinding:
      position: 4

  - id: "#insai"
    type: File
    inputBinding:
      position: 5

  - id: "#infq"
    type: File
    inputBinding:
      position: 6

  - id: "#outsam"
    type: string
    inputBinding:
      position: 1
      prefix: "-f"

  - id: "#max_occ"
    type: ["null", int ]
    description: |
      Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]
    inputBinding:
      position: 1
      prefix: "-n"

  - id: "#RG_line"
    type: ["null", string ]
    description: |
      Specify the read group in a format like `@RG\tID:foo\tSM:bar'. [null]
    inputBinding:
      position: 1
      prefix: "-r"

outputs:
  - id: "#aligned"
    type: File
    outputBinding:
      glob:
        engine: cwl:JsonPointer
        script: /job/outsam

baseCommand: ["bwa","samse"]


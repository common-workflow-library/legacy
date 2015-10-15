#!/usr/bin/env cwl-runner

class: Workflow
requirements:
  - import: ../../tools/envvar-global.cwl

inputs:
  - id: "#input"
    type: File

  - id: "#output"
    type: string

  - id: "#genome1"
    type: string

  - id: "#genome2"
    type: string

outputs:
  - id: "#outfile"
    type: File
    source: "#liftover.output"

steps:
  - id: "#chaindownload"
    run:
      class: CommandLineTool
      requirements:
        - import: node-engine.cwl
      inputs:
        - id: "#genome1"
          type: string

        - id: "#genome2"
          type: string

        - id: "#downloadtemplate"
          type: string
          default: "http://hgdownload.cse.ucsc.edu/goldenPath/{genome1}/liftOver/{genome1}To{genome2}.over.chain.gz"
          inputBinding:
            position: 2
            prefix: "-O"
            valueFrom:
              engine: node-engine.cwl
              script: |
                {
                  var g1=$job['genome1'].toLowerCase();
                  var g2=$job['genome2'].charAt(0).toUpperCase()+$job['genome2'].toLowerCase().slice(1);
                  return $self.replace(/{genome1}/g,g1).replace(/{genome2}/g,g2);
                }
      outputs:
        - id: "#overchain"
          type: File
          outputBinding:
            glob:
              engine: node-engine.cwl
              script: |
                {
                  var g1=$job['genome1'].toLowerCase();
                  var g2=$job['genome2'].charAt(0).toUpperCase()+$job['genome2'].toLowerCase().slice(1);
                  return $job['downloadtemplate'].split('/').slice(-1)[0].replace(/{genome1}/g,g1).replace(/{genome2}/g,g2);
                }

      baseCommand: ["curl","-s"]

    inputs:
      - {id: "#chaindownload.genome1", source: "#genome1"}
      - {id: "#chaindownload.genome2", source: "#genome2"}
    outputs:
      - {id: "#chaindownload.overchain"}

#  - id: "#sort"
#    run: {import: ../../tools/linux-sort.cwl}
#    inputs:
#      - {id: "#sort.input", source: "#genomecov.genomecoverage" }
#      - {id: "#sort.key", default: ["1,1","2,2n"] }
#    outputs:
#      - {id: "#sort.sorted"}

  - id: "liftover"
    run: {import: ../../tools/ucsc-liftOver.cwl}
    inputs:
      - {id: "#liftover.oldFile", source: "#input"}
      - {id: "#liftover.mapChain", source: "#chaindownload.overchain"}
      - {id: "#liftover.newFile", source: "#output" }
      - {id: "#liftover.unMapped", default: 'unMapped.bed' }
    outputs:
      - {id: "#liftover.output"}

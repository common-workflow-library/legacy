#!/usr/bin/env cwl-runner

- id: "#echocmd"
  class: CommandLineTool
  inputs: []
  outputs:
    - id: "#echo-out"
      type: File
      label: "Printed Message"
      description: "The file containing the message"
      outputBinding:
        glob: messageout.txt
  baseCommand: echo
  arguments:
   - Hello
   - world
   - "!"

#!/usr/bin/env cwl-runner

- id: "#echocmd"
  class: CommandLineTool
  inputs:
    - id: "#echo-in"
      type: string
      label: "Message"
      description: "The message to print"
      default: "Hello World"
  outputs:
    - id: "#echo-out"
      type: File
      label: "Printed Message"
      description: "The file containing the message"
      outputBinding:
        glob: messageout.txt
  baseCommand: echo
  stdout: messageout.txt

- id: "#main"
  class: Workflow
  label: "Hello World"
  description: "Puts a message into a file using echo"
  inputs: []
  outputs:
    - id: "#main.output"
      type: File
      source: "#param.output"
  steps :
    - id: "#step0"
      run: {import: "#echocmd"}
      inputs: []
      outputs:
        - { id: "#main.output" }



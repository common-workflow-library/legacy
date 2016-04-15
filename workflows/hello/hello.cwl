#!/usr/bin/env cwl-runner


- id: echocmd
  class: CommandLineTool
  cwlVersion: "cwl:draft-3"
  inputs:
    - id: echo-in
      type: string
      label: "Message"
      description: "The message to print"
      default: "Hello World"
      inputBinding: {}
  outputs:
    - id: echo-out
      type: File
      label: "Printed Message"
      description: "The file containing the message"
      outputBinding:
        glob: messageout.txt
  baseCommand: echo
  stdout: messageout.txt

- id: main
  cwlVersion: "cwl:draft-3"
  class: Workflow
  label: "Hello World"
  description: "Puts a message into a file using echo"
  inputs: []
  outputs:
    - id: output
      type: File
      source: "#main/step0/echo-out"
  steps :
    - id: step0
      run: "#echocmd"
      inputs: []
      outputs:
        - { id: echo-out }

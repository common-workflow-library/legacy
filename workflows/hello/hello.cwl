#!/usr/bin/env cwl-runner

cwlVersion: v1.0
$graph:
- id: echocmd
  class: CommandLineTool
  inputs:
    echo-in:
      type: string
      label: "Message"
      doc: "The message to print"
      default: "Hello World"
      inputBinding: {}
  outputs:
    echo-out:
      type: stdout
      label: "Printed Message"
      doc: "The file containing the message"
  baseCommand: echo
  stdout: messageout.txt

- id: main
  class: Workflow
  label: "Hello World"
  doc: "Puts a message into a file using echo"
  inputs: []
  outputs:
    output:
      type: File
      outputSource: step0/echo-out
  steps:
    step0:
      run: "#echocmd"
      in: []
      out: [echo-out]

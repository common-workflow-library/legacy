#!/usr/bin/env cwl-runner

cwlVersion: v1.0
$graph:
- id: echocmd
  class: CommandLineTool
  inputs:
    echo-in-message:
      type: string
      label: "Message"
      doc: "The message to print"
      default: "Hello World"
      inputBinding: {}
    echo-in-outputfile:
      type: string
      label: "Output file"
      doc: "The file containing the message"
  outputs:
    echo-out-filename:
      type: stdout
      label: "Printed Message"
      doc: "The file containing the message"
  baseCommand: echo
  arguments:
     - "-n"
     - "-e"
  stdout: $(inputs['echo-in-outputfile'])

- id: main
  class: Workflow
  label: "Hello World"
  doc: "Puts a message into a file using echo"
  inputs:
     usermessage: string
     useroutput: string
  outputs:
    output:
      type: File
      outputSource: step0/echo-out-filename
  steps :
    - id: step0
      run: "#echocmd"
      in:
        echo-in-message: usermessage
        echo-in-outputfile: useroutput
      out: [echo-out-filename]

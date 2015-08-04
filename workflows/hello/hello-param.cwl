#!/usr/bin/env cwl-runner


- id: "#echocmd"
  class: CommandLineTool
  inputs:
    - id: "#echo-in-message"
      type: string
      label: "Message"
      description: "The message to print"
      default: "Hello World"
      inputBinding: {}
    - id: "#echo-in-outputfile"
      type: string
      label: "Output file"
      description: "The file containing the message"
  outputs:
    - id: "#echo-out-filename"
      type: File
      label: "Printed Message"
      description: "The file containing the message"
      outputBinding:
          glob:
             engine: cwl:JsonPointer
             script: /job/echo-in-outputfile
  baseCommand: echo
  arguments:
     - "-n"
     - "-e"
  stdout:
    engine: cwl:JsonPointer
    script: /job/echo-in-outputfile


- id: "#main"
  class: Workflow
  label: "Hello World"
  description: "Puts a message into a file using echo"
  inputs:
     - id: "#usermessage"
       type: string
     - id: "#useroutput"
       type: string
  outputs:
    - id: "#main.output"
      type: File
      source: "#echocmd.echo-out-filename"
  steps :
    - id: "#step0"
      run: {import: "#echocmd"}
      inputs:
        - { id: "#echocmd.echo-in-message" ,  source: "#usermessage"  }
        - { id: "#echocmd.echo-in-outputfile",  source: "#useroutput"  }
      outputs:
        - { id: "#echocmd.echo-out-filename", defaul: "default-output.txt" }



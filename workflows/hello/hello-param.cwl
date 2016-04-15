#!/usr/bin/env cwl-runner

- id: echocmd
  class: CommandLineTool
  cwlVersion: cwl:draft-3
  inputs:
    - id: echo-in-message
      type: string
      label: "Message"
      description: "The message to print"
      default: "Hello World"
      inputBinding: {}
    - id: echo-in-outputfile
      type: string
      label: "Output file"
      description: "The file containing the message"
  outputs:
    - id: echo-out-filename
      type: File
      label: "Printed Message"
      description: "The file containing the message"
      outputBinding:
          glob: $(inputs['echo-in-outputfile'])
  baseCommand: echo
  arguments:
     - "-n"
     - "-e"
  stdout: $(inputs['echo-in-outputfile'])


- id: main
  class: Workflow
  label: "Hello World"
  description: "Puts a message into a file using echo"
  inputs:
     - id: usermessage
       type: string
     - id: useroutput
       type: string
  outputs:
    - id: output
      type: File
      source: "#main/step0/echo-out-filename"
  steps :
    - id: step0
      run: "#echocmd"
      inputs:
        - { id: echo-in-message,  source: "#main/usermessage"  }
        - { id: echo-in-outputfile,  source: "#main/useroutput"  }
      outputs:
        - { id: echo-out-filename }


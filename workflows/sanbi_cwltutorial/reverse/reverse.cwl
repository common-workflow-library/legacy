#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: pvanheus/reverse:latest

baseCommand: reverse.py

inputs:
  dnafile:
    type: File
    inputBinding:
      position: 1

stdout: $(inputs.dnafile.nameroot)_reversed$(inputs.dnafile.nameext)

outputs:
  rev_dnafile:
    type: stdout

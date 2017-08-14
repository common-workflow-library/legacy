#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: pvanheus/complement:latest

baseCommand: complement.py

inputs:
  dnafile:
    type: File
    inputBinding:
      position: 1

stdout: $(inputs.dnafile.nameroot)_complement$(inputs.dnafile.nameext)

outputs:
  comp_dnafile:
    type: stdout

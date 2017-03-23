cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: pvanheus/reverse:latest

baseCommand: reverse.py

inputs:
  infile:
    type: File
    inputBinding:
      position: 1

stdout: $(inputs.infile.nameroot)_reversed$(inputs.infile.nameext)

outputs:
  outfile:
    type: stdout

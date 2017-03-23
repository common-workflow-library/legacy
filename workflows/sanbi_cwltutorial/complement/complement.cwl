cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: pvanheus/complement:latest

baseCommand: complement.py

inputs:
  infile:
    type: File
    inputBinding:
      position: 1

stdout: $(inputs.infile.nameroot)_complement$(inputs.infile.nameext)

outputs:
  outfile:
    type: stdout

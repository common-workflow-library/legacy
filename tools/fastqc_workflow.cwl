cwlVersion: v1.0
class: Workflow

inputs:
 reads1: File

outputs:
  zippedFile:
    type: File
    outputSource: fastq/zippedFile
  report:
    type: Directory
    outputSource: fastq/report

steps:
  fastq:
    run: fastqc.cwl
    in:
      fastqFile: reads1
    out: [zippedFile,report]

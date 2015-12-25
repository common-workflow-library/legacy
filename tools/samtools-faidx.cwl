#!/usr/bin/env cwl-runner
#
# To use it as stand alone tool current working directory should not have input .fa file
#    example: "./samtools-faidx.cwl --input=./test-files/mm10.fa"
# As part of a workflow should be no problem at all

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

requirements:
- $import: envvar-global.cwl
- $import: samtools-docker.cwl
- class: InlineJavascriptRequirement
- class: CreateFileRequirement
  fileDef:
  - filename: $(inputs.input.path.split('/').slice(-1)[0])
    fileContent: $(inputs.input)

inputs:
- id: '#input'
  type: File
  description: '<file.fa|file.fa.gz>'

- id: '#region'
  type: ["null",string]
  inputBinding:
    position: 2

outputs:
- id: "#index"
  type: File
  outputBinding:
    glob: $(inputs.input.path.split('/').slice(-1)[0]) #+'.fai')
    secondaryFiles:
    - .fai
    - .gzi

baseCommand:
- samtools
- faidx

arguments:
- valueFrom: $(inputs.input.path.split('/').slice(-1)[0])
  position: 1

description: 'Usage:   samtools faidx <file.fa|file.fa.gz> [<reg> [...]]'

$namespaces:
  schema: http://schema.org/

$schemas:
- https://sparql-test.commonwl.org/schema.rdf

schema:mainEntity:
  $import: samtools-ontology.inc

description: |
  samtools-faidx.cwl is developed for CWL consortium

schema:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-faidx.cwl
schema:codeRepository: https://github.com/common-workflow-language/workflows
schema:license: http://www.apache.org/licenses/LICENSE-2.0

schema:isPartOf:
  class: schema:CreativeWork
  schema:name: "Common Workflow Language"
  schema:url: http://commonwl.org/

schema:author:
  $import: https://scidap.com/porter.yaml

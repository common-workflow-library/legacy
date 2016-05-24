#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

description: |
  samtools-rmdup.cwl is developed for CWL consortium

requirements:
- $import: envvar-global.yml
- $import: samtools-docker.yml
- class: InlineJavascriptRequirement

inputs:
- id: "input"
  type: File
  description: |
    Input bam file.
  inputBinding:
    position: 2

- id: "output_name"
  type: string
  inputBinding:
    position: 3

- id: "single_end"
  type: boolean
  default: false
  description: |
    rmdup for SE reads

- id: "pairend_as_se"
  type: boolean
  default: false
  description: |
    treat PE reads as SE in rmdup (force -s)

outputs:
- id: "rmdup"
  type: File
  description: "File with removed duplicates"
  outputBinding:
    glob: $(inputs.output_name)

baseCommand: ["samtools", "rmdup"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: samtools-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-rmdup.cwl
s:codeRepository: https://github.com/common-workflow-language/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: "Common Workflow Language"
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: "Andrey Kartashov"
  s:email: mailto:Andrey.Kartashov@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  s:worksFor:
  - class: s:Organization
    s:name: "Cincinnati Children's Hospital Medical Center"
    s:location: "3333 Burnet Ave, Cincinnati, OH 45229-3026"
    s:department:
    - class: s:Organization
      s:name: "Barski Lab"

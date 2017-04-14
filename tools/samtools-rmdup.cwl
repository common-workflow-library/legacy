#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- $import: samtools-docker.yml
- class: InlineJavascriptRequirement

inputs:
  single_end:
    type: boolean
    default: false
    doc: |
      rmdup for SE reads
  input:
    type: File
    inputBinding:
      position: 2

    doc: |
      Input bam file.
  output_name:
    type: string
    inputBinding:
      position: 3

  pairend_as_se:
    type: boolean
    default: false
    doc: |
      treat PE reads as SE in rmdup (force -s)
outputs:
  rmdup:
    type: File
    outputBinding:
      glob: $(inputs.output_name)

    doc: File with removed duplicates
baseCommand: [samtools, rmdup]
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
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: Andrey Kartashov
  s:email: mailto:Andrey.Kartashov@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab
doc: |
  samtools-rmdup.cwl is developed for CWL consortium


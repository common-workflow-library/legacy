#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: alea-docker.yml
- class: InlineJavascriptRequirement

inputs:
  hapsDir:
    type: File
    inputBinding:
      position: 2

    doc: |
      path to the directory containing the .haps files
  unphased:
    type: File
    inputBinding:
      position: 3

    doc: |
      path to the vcf file containing unphased SNPs and Indels
  outputPrefix:
    type: string
    inputBinding:
      position: 3

    doc: |
      output file prefix including the path but not the extension
outputs:
  phasevcf:
    type: File
    outputBinding:
      glob: $(inputs.outputPrefix+".vcf.gz")

    doc: Creates the file outputPrefix.vcf.gz
baseCommand: [alea, phaseVCF]
$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: alea-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/alea-phaseVCF.cwl
s:codeRepository: https://github.com/common-workflow-language/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

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


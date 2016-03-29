#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

requirements:
  - $import: envvar-global.cwl
  - $import: alea-docker.cwl
  - class: InlineJavascriptRequirement

inputs:
  - id: "#hapsDir"
    type: File
    description: |
      path to the directory containing the .haps files
    inputBinding:
      position: 2

  - id: "#unphased"
    type: File
    description: |
      path to the vcf file containing unphased SNPs and Indels
    inputBinding:
      position: 3

  - id: "#outputPrefix"
    type: string
    description: |
      output file prefix including the path but not the extension
    inputBinding:
      position: 3

outputs:
  - id: "#phasevcf"
    type: File
    description: "Creates the file outputPrefix.vcf.gz"
    outputBinding:
      glob: $(inputs.outputPrefix+".vcf.gz")

baseCommand: ["alea", "phaseVCF"]

$namespaces:
  schema: http://schema.org/

$schemas:
- https://sparql-test.commonwl.org/schema.rdf

schema:mainEntity:
  $import: alea-metadata.yaml

schema:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/alea-phaseVCF.cwl
schema:codeRepository: https://github.com/common-workflow-language/workflows
schema:license: http://www.apache.org/licenses/LICENSE-2.0
schema:isPartOf:
  class: schema:CreativeWork
  schema:name: "Common Workflow Language"
  schema:url: http://commonwl.org/

schema:author:
  class: schema:Person
  schema:name: "Andrey Kartashov"
  schema:email: mailto:Andrey.Kartashov@cchmc.org
  schema:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  schema:worksFor:
  - class: schema:Organization
    schema:name: "Cincinnati Children's Hospital Medical Center"
    schema:location: "3333 Burnet Ave, Cincinnati, OH 45229-3026"
    schema:department:
    - class: schema:Organization
      schema:name: "Barski Lab"

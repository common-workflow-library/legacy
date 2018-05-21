#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: alea-docker.yml
- class: InlineJavascriptRequirement

inputs:
  strain:
    type: string
    inputBinding:
      separate: false
      prefix: --strain=
      position: 4

    doc: name of strain1 exactly as specified in the vcf file (e.g. hap1)
  phased:
    type: File
    inputBinding:
      separate: false
      prefix: --input-vcf=
      position: 3
    secondaryFiles:
    - .tbi
    doc: |
      the phased variants vcf file (SNPs or/and Indels)
  reference:
    type: File
    inputBinding:
      separate: false
      prefix: --input-fasta=
      position: 2
    secondaryFiles:
    - .fai
    doc: |
      the reference genome fasta file
  filename:
    type: string
    inputBinding:
      position: 5
      separate: false
      prefix: --output-fasta=
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
    secondaryFiles:
    - .fai
    - .refmap
baseCommand: [java, -Xms4G, -Xmx8G, -jar, /usr/local/bin/alea.jar, insilico]
$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html
#- http://topbraid.org/schema/schema.rdf

s:mainEntity:
  $import: alea-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/alea-insilico.cwl
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


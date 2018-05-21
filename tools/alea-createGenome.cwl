#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: alea-docker.yml
- class: EnvVarRequirement
  envDef:
  - envName: AL_USE_CONCATENATED_GENOME
    envValue: $(inputs.CONCATENATED_GENOME?"1":"0")
  - envName: AL_BWA_ALN_PARAMS
    envValue: -k 0 -n 0 -t 4
  - envName: AL_DIR_TOOLS
    envValue: /usr/local/bin/
- class: InlineJavascriptRequirement

inputs:
  CONCATENATED_GENOME:
    type: boolean
    default: false

  reference:
    type: File
    inputBinding:
      position: 2
    secondaryFiles:
    - .fai
    doc: |
      the reference genome fasta file
  phasedindels:
    type: File?
    inputBinding:
      position: 4
    secondaryFiles:
    - .tbi
    doc: |
      the phased Indels (should be specified second)
  strain2:
    type: string
    inputBinding:
      position: 6

    doc: |
      name of strain2 exactly as specified in the vcf file (e.g. hap2)
  strain1:
    type: string
    inputBinding:
      position: 5

    doc: name of strain1 exactly as specified in the vcf file (e.g. hap1)
  phased:
    type: File
    inputBinding:
      position: 3
    secondaryFiles:
    - .tbi
    doc: |
      the phased variants vcf file (including SNPs and Indels)
      or the phased SNPs (should be specified first)
  outputDir:
    type: string
    inputBinding:
      position: 7

    doc: |
      location of the output directory
outputs:
  strain12_indices:
    type: File?
    outputBinding:
      glob: $(inputs.CONCATENATED_GENOME?inputs.outputDir+"/"+inputs.strain1+"_"+inputs.strain2+".fasta":[])
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .sa
  strain1_indices:
    type: File
    outputBinding:
      glob: $(inputs.outputDir+"/"+inputs.strain1+".fasta")
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .refmap
    - .sa
  strain2_indices:
    type: File
    outputBinding:
      glob: $(inputs.outputDir+"/"+inputs.strain1+".fasta")
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .refmap
    - .sa
baseCommand: [alea, createGenome]
arguments:
- valueFrom: $(inputs.phasedindels?"-snps-indels-separately":[])
  position: 1 

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: alea-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/alea-createGenome.cwl
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


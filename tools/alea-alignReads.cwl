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

  outputPrefix:
    type: string
    inputBinding:
      position: 7

    doc: |
      location of the output directory
  input_reads:
    type:
      type: array
      items: File
    inputBinding:
      position: 2

    doc: |
      input_reads_1 the 1st input reads file in fastq.
      (fastq.gz or bam is supported when using BWA)
      input_reads_2 (paired end) the 2nd input reads file in fastq. (fastq.gz or bam is supported when using BWA)
  genome2:
    type: File?
    inputBinding:
      position: 3
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .pac
    - .sa
    doc: |
      (when AL_USE_CONCATENATED_GENOME=0)
      path to the indexed reference for 2nd insilico genome (of strain2). for BWA, specifiy the fasta file.
      for Bowtie, specify basename of index files.
  genome1:
    type: File
    inputBinding:
      position: 3
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .pac
    - .sa
    doc: |
      (when AL_USE_CONCATENATED_GENOME=0)
      path to the indexed reference for 1st insilico genome (of strain1).
      for BWA, specifiy the fasta file.
      for Bowtie, specify index filename prefix (minus trailing .X.ebwt or .X.bt2)
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
outputs: {}
#- id: "strain1_indices"
#  type: File
#  outputBinding:
#    glob: $(inputs.outputDir+"/"+inputs.strain1+".fasta")
#    secondaryFiles:
#    - ".amb"
#    - ".ann"
#    - ".bwt"
#    - ".fai"
#    - ".pac"
#    - ".refmap"
#    - ".sa"

baseCommand: [alea, alignReads]
arguments:
- valueFrom: $(inputs.input_reads.length==1?"-s":"-p")
  position: 1 

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: alea-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/alea-alignReads.cwl
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


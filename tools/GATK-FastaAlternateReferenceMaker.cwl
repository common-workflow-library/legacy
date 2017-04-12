#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- $import: GATK-docker.yml
- class: InlineJavascriptRequirement

inputs:
  output_filename:
    type: string
    inputBinding:
      prefix: -o
      position: 4

  snpmask:
    type: File?
    inputBinding:
      prefix: --snpmask
      position: 4

    doc: |
      SNP mask VCF file
  intervals:
    type: File?
    inputBinding:
      prefix: -L
      position: 4

  vcf:
    type: File
    inputBinding:
      prefix: -V
      position: 4

    doc: |
      Input VCF file
      Variants from this VCF file are used by this tool as input. The file must at least contain the standard VCF header lines, but can be empty (i.e., no variants are contained in the file).
      --variant binds reference ordered data. This argument supports ROD files of the following types: BCF2, VCF, VCF3
  reference:
    type: File
    inputBinding:
      prefix: -R
      position: 4

    doc: |
      Input reference fasta or fasta.gz
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

    doc: 'An output file created by the walker. Will overwrite contents if file exists

      '
baseCommand: [java]
arguments:
- valueFrom: -Xmx4g
  position: 1
- valueFrom: /usr/local/bin/GenomeAnalysisTK.jar
  position: 2
  prefix: -jar
- valueFrom: FastaAlternateReferenceMaker
  position: 3
  prefix: -T
$namespaces:
  schema: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

schema:mainEntity:
  $import: GATK-metadata.yaml

schema:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-FastaAlternateReferenceMaker.cwl
schema:codeRepository: https://github.com/common-workflow-language/workflows
schema:license: http://www.apache.org/licenses/LICENSE-2.0
schema:isPartOf:
  class: schema:CreativeWork
  schema:name: Common Workflow Language
  schema:url: http://commonwl.org/

schema:author:
  class: schema:Person
  schema:name: Andrey Kartashov
  schema:email: mailto:Andrey.Kartashov@cchmc.org
  schema:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  schema:worksFor:
  - class: schema:Organization
    schema:name: Cincinnati Children's Hospital Medical Center
    schema:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    schema:department:
    - class: schema:Organization
      schema:name: Barski Lab
doc: |
  GATK-FastaAlternateReferenceMaker.cwl is developed for CWL consortium


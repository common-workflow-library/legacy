#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- $import: samtools-docker.yml
- class: InlineJavascriptRequirement

inputs:
  compression_level:
    type: int?
    inputBinding:
      prefix: -l
    doc: |
      Set compression level, from 0 (uncompressed) to 9 (best)
  threads:
    type: int?
    inputBinding:
      prefix: -@

    doc: Set number of sorting and compression threads [1]
  memory:
    type: string?
    inputBinding:
      prefix: -m
    doc: |
      Set maximum memory per thread; suffix K/M/G recognized [768M]
  input:
    type: File
    inputBinding:
      position: 1

    doc: Input bam file.
  output_name:
    type: string
    inputBinding:
      position: 2

    doc: Desired output filename.
  sort_by_name:
    type: boolean?
    inputBinding:
      prefix: -n

    doc: Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
outputs:
  sorted:
    type: File
    outputBinding:
      glob: $(inputs.output_name)

baseCommand: [samtools, sort]
arguments:
- -f
$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: samtools-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-sort.cwl
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
  samtools-sort.cwl is developed for CWL consortium
    Usage: samtools sort [options...] [in.bam]
    Options:
      -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
      -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
      -n         Sort by read name
      -o FILE    Write final output to FILE rather than standard output
      -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')   (either -O or
      -T PREFIX  Write temporary files to PREFIX.nnnn.bam       -T is required)
      -@ INT     Set number of sorting and compression threads [1]

    Legacy usage: samtools sort [options...] <in.bam> <out.prefix>
    Options:
      -f         Use <out.prefix> as full final filename rather than prefix
      -o         Write final output to stdout rather than <out.prefix>.bam
      -l,m,n,@   Similar to corresponding options above



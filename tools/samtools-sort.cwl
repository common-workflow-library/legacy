#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

description: |
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

requirements:
- $import: envvar-global.yml
- $import: samtools-docker.yml
- class: InlineJavascriptRequirement

inputs:
- id: "compression_level"
  type: ["null", int]
  description: |
    Set compression level, from 0 (uncompressed) to 9 (best)
  inputBinding:
    prefix: "-l"

- id: "memory"
  type: ["null", string]
  description: |
    Set maximum memory per thread; suffix K/M/G recognized [768M]
  inputBinding:
    prefix: "-m"

- id: "sort_by_name"
  type: ["null", boolean]
  description: "Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates."
  inputBinding:
    prefix: -n

- id: "threads"
  type: ["null", int]
  description: "Set number of sorting and compression threads [1]"
  inputBinding:
    prefix: -@

- id: "output_name"
  type: string
  description: "Desired output filename."
  inputBinding:
    position: 2

- id: "input"
  type: File
  description:
    Input bam file.
  inputBinding:
    position: 1

outputs:
- id: "sorted"
  type: File
  outputBinding:
    glob: $(inputs.output_name)

baseCommand: ["samtools", "sort"]

arguments:
  - "-f"

$namespaces:
  s: http://schema.org/

$schemas:
- https://sparql-test.commonwl.org/schema.rdf

s:mainEntity:
  $import: samtools-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-sort.cwl
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

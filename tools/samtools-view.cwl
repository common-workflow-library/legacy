#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- $import: samtools-docker.yml
- class: InlineJavascriptRequirement

inputs:
  isbam:
    type: boolean
    default: false
    inputBinding:
      position: 2
      prefix: -b
    doc: |
      output in BAM format
  readswithoutbits:
    type: int?
    inputBinding:
      position: 1
      prefix: -F
    doc: |
      only include reads with none of the bits set in INT set in FLAG [0]
  collapsecigar:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -B
    doc: |
      collapse the backward CIGAR operation
  readsingroup:
    type: string?
    inputBinding:
      position: 1
      prefix: -r
    doc: |
      only include reads in read group STR [null]
  bedoverlap:
    type: File?
    inputBinding:
      position: 1
      prefix: -L
    doc: |
      only include reads overlapping this BED FILE [null]
  uncompressed:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -u
    doc: |
      uncompressed BAM output (implies -b)
  readtagtostrip:
    type: string[]?
    inputBinding:
      position: 1

    doc: |
      read tag to strip (repeatable) [null]
  input:
    type: File
    inputBinding:
      position: 4

    doc: |
      Input bam file.
  readsquality:
    type: int?
    inputBinding:
      position: 1
      prefix: -q
    doc: |
      only include reads with mapping quality >= INT [0]
  readswithbits:
    type: int?
    inputBinding:
      position: 1
      prefix: -f
    doc: |
      only include reads with all bits set in INT set in FLAG [0]
  cigar:
    type: int?
    inputBinding:
      position: 1
      prefix: -m
    doc: |
      only include reads with number of CIGAR operations
      consuming query sequence >= INT [0]
  iscram:
    type: boolean
    default: false
    inputBinding:
      position: 2
      prefix: -C
    doc: |
      output in CRAM format
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: -@
    doc: |
      number of BAM compression threads [0]
  fastcompression:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: '-1'
    doc: |
      use fast BAM compression (implies -b)
  samheader:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -h
    doc: |
      include header in SAM output
  count:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -c
    doc: |
      print only the count of matching records
  randomseed:
    type: float?
    inputBinding:
      position: 1
      prefix: -s
    doc: |
      integer part sets seed of random number generator [0];
      rest sets fraction of templates to subsample [no subsampling]
  referencefasta:
    type: File?
    inputBinding:
      position: 1
      prefix: -T
    doc: |
      reference sequence FASTA FILE [null]
  region:
    type: string?
    inputBinding:
      position: 5

    doc: |
      [region ...]
  readsingroupfile:
    type: File?
    inputBinding:
      position: 1
      prefix: -R
    doc: |
      only include reads with read group listed in FILE [null]
  readsinlibrary:
    type: string?
    inputBinding:
      position: 1
      prefix: -l
    doc: |
      only include reads in library STR [null]
  output_name:
    type: string
    inputBinding:
      position: 2
      prefix: -o
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_name)

baseCommand: [samtools, view]
$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: samtools-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-view.cwl
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
  samtools-view.cwl is developed for CWL consortium
    Usage:   samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]

    Options: -b       output BAM
             -C       output CRAM (requires -T)
             -1       use fast BAM compression (implies -b)
             -u       uncompressed BAM output (implies -b)
             -h       include header in SAM output
             -H       print SAM header only (no alignments)
             -c       print only the count of matching records
             -o FILE  output file name [stdout]
             -U FILE  output reads not selected by filters to FILE [null]
             -t FILE  FILE listing reference names and lengths (see long help) [null]
             -T FILE  reference sequence FASTA FILE [null]
             -L FILE  only include reads overlapping this BED FILE [null]
             -r STR   only include reads in read group STR [null]
             -R FILE  only include reads with read group listed in FILE [null]
             -q INT   only include reads with mapping quality >= INT [0]
             -l STR   only include reads in library STR [null]
             -m INT   only include reads with number of CIGAR operations
                      consuming query sequence >= INT [0]
             -f INT   only include reads with all bits set in INT set in FLAG [0]
             -F INT   only include reads with none of the bits set in INT
                      set in FLAG [0]
             -x STR   read tag to strip (repeatable) [null]
             -B       collapse the backward CIGAR operation
             -s FLOAT integer part sets seed of random number generator [0];
                      rest sets fraction of templates to subsample [no subsampling]
             -@ INT   number of BAM compression threads [0]

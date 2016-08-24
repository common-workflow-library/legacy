#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

description: |
  ucsc-liftOver.cwl is developed for CWL consortium
    usage:
       liftOver oldFile map.chain newFile unMapped
    oldFile and newFile are in bed format by default, but can be in GFF and
    maybe eventually others with the appropriate flags below.
    The map.chain file has the old genome as the target and the new genome
    as the query.

    ***********************************************************************
    WARNING: liftOver was only designed to work between different
             assemblies of the same organism. It may not do what you want
             if you are lifting between different organisms. If there has
             been a rearrangement in one of the species, the size of the
             region being mapped may change dramatically after mapping.
    ***********************************************************************

requirements:
  - class: InlineJavascriptRequirement
  - $import: envvar-global.yml
  - $import: ucsc-userapps-docker.yml

inputs:
  - id: "#oldFile"
    type: File
    inputBinding:
      position: 2

  - id: "#mapChain"
    type: File
    description: |
      The map.chain file has the old genome as the target and the new genome
      as the query.
    inputBinding:
      position: 3

  - id: "#newFile"
    type: string
    inputBinding:
      position: 4

  - id: "#unMapped"
    type: string
    inputBinding:
      position: 5

  - id: "#gff"
    type: ["null",boolean]
    description: |
      File is in gff/gtf format.  Note that the gff lines are converted
       separately.  It would be good to have a separate check after this
       that the lines that make up a gene model still make a plausible gene
       after liftOver
    inputBinding:
      position: 1
      prefix: "-gff"

  - id: "#genePred"
    type: ["null",boolean]
    description: |
      File is in genePred format
    inputBinding:
      position: 1
      prefix: "-genePred"

  - id: "#sample"
    type: ["null",boolean]
    description: |
      File is in sample format
    inputBinding:
      position: 1
      prefix: "-sample"

  - id: "#bedPlus"
    type: ["null",int]
    description: |
      =N - File is bed N+ format
    inputBinding:
      separate: false
      position: 1
      prefix: "-bedPlus="

  - id: "#positions"
    type: ["null",boolean]
    description: |
      File is in browser "position" format
    inputBinding:
      position: 1
      prefix: "-positions"

  - id: "#hasBin"
    type: ["null",boolean]
    description: |
      File has bin value (used only with -bedPlus)
    inputBinding:
      position: 1
      prefix: "-hasBin"

  - id: "#minMatch"
    type: ["null",int]
    description: |
      -minMatch=0.N Minimum ratio of bases that must remap. Default 0.95
    inputBinding:
      separate: false
      position: 1
      prefix: "-minMatch="

  - id: "#tab"
    type: ["null",boolean]
    inputBinding:
      position: 1
      prefix: "-tab"

  - id: "#pslT"
    type: ["null",boolean]
    description: |
      File is in psl format, map target side only
    inputBinding:
      position: 1
      prefix: "-pslT"

  - id: "#ends"
    type: ["null",int]
    description: |
      =N - Lift the first and last N bases of each record and combine the
               result. This is useful for lifting large regions like BAC end pairs.
    inputBinding:
      separate: false
      position: 1
      prefix: "-ends="

  - id: "#minBlocks"
    type: ["null",int]
    description: |
      .N Minimum ratio of alignment blocks or exons that must map
                    (default 1.00)
    inputBinding:
      separate: false
      position: 1
      prefix: "-minBlocks="

  - id: "#fudgeThick"
    type: ["null",boolean]
    description: |
      (bed 12 or 12+ only) If thickStart/thickEnd is not mapped,
                    use the closest mapped base.  Recommended if using
                    -minBlocks.
    inputBinding:
      position: 1
      prefix: "-fudgeThick"

  - id: "#multiple"
    type: ["null",boolean]
    description: |
      Allow multiple output regions
    inputBinding:
      position: 1
      prefix: "-multiple"

  - id: "#minChainT"
    type: ["null",int]
    description: |
      Minimum chain size in target/query, when mapping
                             to multiple output regions (default 0, 0)
    inputBinding:
      position: 1
      prefix: "-minChainT"

  - id: "#minChainQ"
    type: ["null",int]
    description: |
      Minimum chain size in target/query, when mapping
                             to multiple output regions (default 0, 0)
    inputBinding:
      position: 1
      prefix: "-minChainQ"

  - id: "#minSizeQ"
    type: ["null",int]
    description: |
      Min matching region size in query with -multiple.
    inputBinding:
      position: 1
      prefix: "-minSizeQ"

  - id: "#chainTable"
    type: ["null",string]
    description: |
      Min matching region size in query with -multiple.
    inputBinding:
      position: 1
      prefix: "-chainTable"

outputs:
  - id: "#output"
    type: File
    description: "The sorted file"
    outputBinding:
      glob: $(inputs.newFile)

  - id: "#unMappedFile"
    type: File
    description: "The sorted file"
    outputBinding:
      glob: $(inputs.unMapped)

baseCommand: ["liftOver"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ucsc-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/ucsc-liftOver.cwl
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

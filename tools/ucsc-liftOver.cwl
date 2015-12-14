#!/usr/bin/env cwl-runner

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  adms: http://www.w3.org/ns/adms#
  dcat: http://www.w3.org/ns/dcat#

$schemas:
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

adms:includedAsset:
  doap:name: "UCSC userApps"
  doap:description: |
    UCSC genome browser 'kent' bioinformatic utilities
    These are only the command line bioinformatic utilities
    from the kent source tree.
    liftOver - Move annotations from one assembly to another
  doap:homepage: "https://genome.ucsc.edu/util.html"
  dcat:downloadURL: "http://hgdownload.cse.ucsc.edu/admin/exe/userApps.v325.src.tgz"
  doap:release:
  - class: doap:Version
    doap:revision: "v325"
  doap:license: "GPL"
  doap:category: "commandline tool"
  doap:programming-language: "C"
  foaf:publications:
  - id: urn:pmid:20639541
    foaf:title: "(Kent et al., 2010) BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics."
    foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/20639541"
  doap:developer:
  - class: foaf:Organization
    foaf:name: "CIRM Stem Cell Genomics Data Management Center"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Jim Kent"
      foaf:mbox: "mailto:kent@soe.ucsc.edu"

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

doap:name: "ucsc-liftOver.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/ucsc-liftOver.cwl"

dct:isPartOf:
  doap:name: "CWL Workflows"
  doap:homepage: "http://commonwl.org/"
  doap:license: "Apache2"

  doap:implements:
  - class: doap:Specification
    doap:homepage: "http://common-workflow-language.github.io/draft-3/"

  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/common-workflow-language/workflows"

  dct:creator:
  - class: foaf:Organization
    foaf:name: "Curoverse"
    foaf:member:
    - class: foaf:Person
      id: "http://orcid.org/0000-0003-3566-7705"
      foaf:name: "Peter Amstutz"
      foaf:mbox: "mailto:peter.amstutz@curoverse.com"
  - class: foaf:Organization
    foaf:name: "Seven Bridges Genomics"
    foaf:member:
    - class: foaf:Person
      id: "mailto:nebojsa.tijanic@sbgenomics.com"
      foaf:name: "Nebojša Tijanić"
      foaf:mbox: "mailto:nebojsa.tijanic@sbgenomics.com"

  dct:contributor:
  - class: foaf:Organization
    foaf:name: "Seven Bridges Genomics"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Luka Stojanovic"
      foaf:mbox: "mailto:luka.stojanovic@sbgenomics.com"
  - class: foaf:Organization
    foaf:name: "Galaxy Project, Pennsylvania State University"
    foaf:member:
    - class: foaf:Person
      foaf:name: "John Chilton"
      foaf:mbox: "mailto:jmchilton@gmail.com"
  - class: foaf:Organization
    foaf:name: "University of California, Davis"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Michael R. Crusoe"
      foaf:mbox: "mailto:crusoe@ucdavis.edu"
  - class: foaf:Organization
    foaf:name: "Institut Pasteur"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Hervé Ménager"
      foaf:mbox: "mailto:herve.menager@gmail.com"
  - class: foaf:Organization
    foaf:name: "BioDatomics"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Maxim Mikheev"
      foaf:mbox: "mailto:mikhmv@biodatomics.com"
  - class: foaf:Organization
    foaf:name: "University of Manchester"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Stian Soiland-Reyes"
      foaf:mbox: "mailto:soiland-reyes@cs.manchester.ac.uk"

doap:maintainer:
- class: foaf:Organization
  foaf:name: "Barski Lab, Cincinnati Children's Hospital Medical Center"
  foaf:member:
  - class: foaf:Person
    id: "http://orcid.org/0000-0001-9102-5681"
    foaf:openid: "http://orcid.org/0000-0001-9102-5681"
    foaf:name: "Andrey Kartashov"
    foaf:mbox: "mailto:Andrey.Kartashov@cchmc.org"

requirements:
  - class: InlineJavascriptRequirement
  - $import: envvar-global.cwl
  - $import: ucsc-userapps-docker.cwl

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



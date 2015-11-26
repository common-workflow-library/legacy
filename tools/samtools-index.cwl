#!/usr/bin/env cwl-runner

"@context":
  "foaf": "http://xmlns.com/foaf/0.1/"
  "doap": "http://usefulinc.com/ns/doap"
  "adms": "http://purl.org/adms/"
  "admssw": "http://purl.org/adms/sw/"

adms:Asset:
  admssw:SoftwareProject:
    doap:name: "samtools"
    doap:description: >
      A suite of programs for interacting with high-throughput sequencing data.
      It consists of three separate repositories: Samtools (Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format),
      BCFtools (Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants)
      and HTSlib (A C library for reading/writing high-throughput sequencing data).
    doap:homepage: "http://www.htslib.org/"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/samtools/samtools.git"
    doap:release:
      - doap:revision: "1.2-216-gdffc67f"
    doap:license: "MIT, BSD License"
    doap:category: "commandline tool"
    doap:programming-language: "C, Perl"
    foaf:publications:
      - foaf:title: "(Li, 2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics."
        foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/21903627"
      - foaf:title: "(Li, 2011) Improving SNP discovery by base alignment quality. Bioinformatics."
        foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/21320865"
      - foaf:title: "(Li et al., 2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics."
        foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/19505943"
    doap:mailing-list:
      - doap:location: "https://lists.sourceforge.net/lists/listinfo/samtools-help"
      - doap:location: "https://lists.sourceforge.net/lists/listinfo/samtools-devel"
    doap:developer:
      foaf:Group:
        - foaf:Person:
          foaf:name: "Heng Li"
          foaf:Organization: "Sanger Institute"
          doap:description: "wrote most of the initial source codes of SAMtools and various converters."
        - foaf:Person:
          foaf:name: "Bob Handsaker"
          foaf:Organization: "Broad Institute"
          doap:description: |
            A major contributor to the
            SAM/BAM specification. He designed and implemented the BGZF format, the
            underlying indexable compression format for the BAM format. BGZF does
            not support arithmetic between file offsets.
        - foaf:Person:
          foaf:name: "Jue Ruan"
          foaf:Organization: "Beijing Genome Institute"
          doap:description: |
            Designed and implemented the
            RAZF format, an alternative indexable compression format. RAZF is no longer
            used by or provided with SAMtools. Source code remains available in older
            SAMtools 0.1.x releases and from the standalone branch in the repository.
        - foaf:Person:
          foaf:name: "Colin Hercus"
          doap:description: "updated novo2sam.pl to support gapped alignment by novoalign."
        - foaf:Person:
          foaf:name: "Petr Danecek"
          doap:description: "contributed the header parsing library sam_header.c and sam2vcf.pl script."
  adms:AssetDistribution:
    doap:name: "samtools-index.cwl"
    doap:description: "Developed for CWL consortium http://commonwl.org/"
    doap:specification: "http://common-workflow-language.github.io/draft-3/"
    doap:release: "cwl:draft-3.dev2"
    doap:homepage: "http://commonwl.org/"
    doap:location: "https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-index.cwl"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/common-workflow-language/workflows"
    doap:maintainer:
      foaf:Person:
        foaf:openid: "http://orcid.org/0000-0001-9102-5681"
        foaf:name: "Andrey Kartashov"
        foaf:mbox: "mailto:Andrey.Kartashov@cchmc.org"
        foaf:organization: "Cincinnati Children's Hospital Medical Center"

cwlVersion: "cwl:draft-3.dev2"

class: CommandLineTool

description: |
  Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
  Options:
    -b       Generate BAI-format index for BAM files [default]
    -c       Generate CSI-format index for BAM files
    -m INT   Set minimum interval size for CSI indices to 2^INT [14]

requirements:
  - "@import": envvar-global.cwl
  - "@import": samtools-docker.cwl
  - class: InlineJavascriptRequirement
    expressionLib:
    - "var new_ext = function() { var ext=inputs.bai?'.bai':inputs.csi?'.csi':'.bai'; return inputs.input.path.split('/').slice(-1)[0]+ext; };"
#  - class: CreateFileRequirement
#    fileDef:
#      - filename: $(inputs.input.path.split('/').slice(-1)[0])
#        fileContent: $(inputs.input)

inputs:
  - id: "#input"
    type: File
    description: |
      Input bam file.
    inputBinding:
      position: 2

  - id: "#bai"
    type: boolean
    default: false
    description: |
      Generate BAI-format index for BAM files [default]

  - id: "#csi"
    type: boolean
    default: false
    description: |
      Generate CSI-format index for BAM files

  - id: "#interval"
    type: ["null", int]
    description: |
      Generate CSI-format index for BAM files
    inputBinding:
      position: 1
      prefix: "-m"

outputs:
  - id: "#sorted"
    type: File
    description: "The sorted file"
    outputBinding:
      glob: $(new_ext())

baseCommand: ["samtools", "index"]

arguments:
  - valueFrom: $(inputs.bai?'-b':inputs.csi?'-c':[])
    position: 1
  - valueFrom: $(new_ext())
    position: 3

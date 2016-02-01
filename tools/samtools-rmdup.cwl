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
  doap:name: "samtools"
  doap:description: >
    A suite of programs for interacting with high-throughput sequencing data.
    It consists of three separate repositories: Samtools (Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format),
    BCFtools (Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants)
    and HTSlib (A C library for reading/writing high-throughput sequencing data).
  doap:homepage: "http://www.htslib.org/"
  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/samtools/samtools.git"
  doap:release:
  - class: doap:Version
    doap:revision: "1.2-242-4d56437"
  doap:license: "MIT, BSD License"
  doap:category: "commandline tool"
  doap:programming-language: "C, Perl"
  foaf:publications:
  - id: urn:pmid:21903627
    foaf:title: "(Li, 2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics."
    foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/21903627"
  - id: urn:pmid:21320865
    foaf:title: "(Li, 2011) Improving SNP discovery by base alignment quality. Bioinformatics."
    foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/21320865"
  - id: urn:pmid:19505943
    foaf:title: "(Li et al., 2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics."
    foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/19505943"
  doap:mailing-list:
  - doap:location: "https://lists.sourceforge.net/lists/listinfo/samtools-help"
  - doap:location: "https://lists.sourceforge.net/lists/listinfo/samtools-devel"
  doap:developer:
  - class: foaf:Organization
    foaf:name: "Sanger Institute"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Heng Li"
      doap:description: "wrote most of the initial source codes of SAMtools and various converters."
  - class: foaf:Organization
    foaf:name: "Broad Institute"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Bob Handsaker"
      doap:description: |
        A major contributor to the
        SAM/BAM specification. He designed and implemented the BGZF format, the
        underlying indexable compression format for the BAM format. BGZF does
        not support arithmetic between file offsets.
  - class: foaf:Organization
    foaf:name: "Beijing Genome Institute"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Jue Ruan"
      doap:description: |
        Designed and implemented the
        RAZF format, an alternative indexable compression format. RAZF is no longer
        used by or provided with SAMtools. Source code remains available in older
        SAMtools 0.1.x releases and from the standalone branch in the repository.
  - class: foaf:Person
    foaf:name: "Colin Hercus"
    doap:description: "updated novo2sam.pl to support gapped alignment by novoalign."
  - class: foaf:Person
    foaf:name: "Petr Danecek"
    doap:description: "contributed the header parsing library sam_header.c and sam2vcf.pl script."

description: |
  samtools-rmdup.cwl is developed for CWL consortium

doap:name: "samtools-rmdup.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-rmdup.cwl"

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
- $import: envvar-global.cwl
- $import: samtools-docker.cwl
- class: InlineJavascriptRequirement

inputs:
- id: "input"
  type: File
  description: |
    Input bam file.
  inputBinding:
    position: 2

- id: "output_name"
  type: string
  inputBinding:
    position: 3

- id: "single_end"
  type: boolean
  default: false
  description: |
    rmdup for SE reads

- id: "pairend_as_se"
  type: boolean
  default: false
  description: |
    treat PE reads as SE in rmdup (force -s)

outputs:
- id: "rmdup"
  type: File
  description: "File with removed duplicates"
  outputBinding:
    glob: $(inputs.output_name)

baseCommand: ["samtools", "rmdup"]

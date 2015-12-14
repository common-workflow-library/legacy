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

doap:name: "samtools-view.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-view.cwl"

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
- $import: envvar-global.cwl
- $import: samtools-docker.cwl
- class: InlineJavascriptRequirement

inputs:
- id: "#input"
  type: File
  description: |
    Input bam file.
  inputBinding:
    position: 4

- id: "#region"
  type: ["null",string]
  description: |
    [region ...]
  inputBinding:
    position: 5

- id: "#output_name"
  type: string
  inputBinding:
    position: 2
    prefix: "-o"

- id: "#isbam"
  type: boolean
  default: false
  description: |
    output in BAM format
  inputBinding:
    position: 2
    prefix: "-b"

- id: "#iscram"
  type: boolean
  default: false
  description: |
    output in CRAM format
  inputBinding:
    position: 2
    prefix: "-C"

- id: "#fastcompression"
  type: boolean
  default: false
  description: |
    use fast BAM compression (implies -b)
  inputBinding:
    position: 1
    prefix: "-1"

- id: "#uncompressed"
  type: boolean
  default: false
  description: |
    uncompressed BAM output (implies -b)
  inputBinding:
    position: 1
    prefix: "-u"

- id: "#samheader"
  type: boolean
  default: false
  description: |
    include header in SAM output
  inputBinding:
    position: 1
    prefix: "-h"

- id: "#count"
  type: boolean
  default: false
  description: |
    print only the count of matching records
  inputBinding:
    position: 1
    prefix: "-c"

- id: "#referencefasta"
  type: ["null",File]
  description: |
    reference sequence FASTA FILE [null]
  inputBinding:
    position: 1
    prefix: "-T"

- id: "#bedoverlap"
  type: ["null",File]
  description: |
    only include reads overlapping this BED FILE [null]
  inputBinding:
    position: 1
    prefix: "-L"

- id: "#readsingroup"
  type: ["null",string]
  description: |
    only include reads in read group STR [null]
  inputBinding:
    position: 1
    prefix: "-r"

- id: "#readsingroupfile"
  type: ["null",File]
  description: |
    only include reads with read group listed in FILE [null]
  inputBinding:
    position: 1
    prefix: "-R"

- id: "#readsquality"
  type: ["null",int]
  description: |
    only include reads with mapping quality >= INT [0]
  inputBinding:
    position: 1
    prefix: "-q"

- id: "#readsinlibrary"
  type: ["null",string]
  description: |
    only include reads in library STR [null]
  inputBinding:
    position: 1
    prefix: "-l"

- id: "#cigar"
  type: ["null",int]
  description: |
    only include reads with number of CIGAR operations
    consuming query sequence >= INT [0]
  inputBinding:
    position: 1
    prefix: "-m"

- id: "#readswithbits"
  type: ["null",int]
  description: |
    only include reads with all bits set in INT set in FLAG [0]
  inputBinding:
    position: 1
    prefix: "-f"

- id: "#readswithoutbits"
  type: ["null",int]
  description: |
    only include reads with none of the bits set in INT set in FLAG [0]
  inputBinding:
    position: 1
    prefix: "-F"

- id: "#readtagtostrip"
  type:
  - "null"
  - type: array
    items: string
    inputBinding:
      prefix: "-x"
  description: |
    read tag to strip (repeatable) [null]
  inputBinding:
    position: 1

- id: "#collapsecigar"
  type: boolean
  default: false
  description: |
    collapse the backward CIGAR operation
  inputBinding:
    position: 1
    prefix: "-B"

- id: "#randomseed"
  type: ["null",float]
  description: |
    integer part sets seed of random number generator [0];
    rest sets fraction of templates to subsample [no subsampling]
  inputBinding:
    position: 1
    prefix: "-s"

- id: "#threads"
  type: ["null",int]
  description: |
    number of BAM compression threads [0]
  inputBinding:
    position: 1
    prefix: "-@"

outputs:
- id: "#output"
  type: File
  outputBinding:
    glob: $(inputs.output_name)

baseCommand: ["samtools", "view"]

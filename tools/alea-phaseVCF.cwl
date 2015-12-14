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
  doap:name: "alea"
  doap:description: >
    ALEA is a computational toolbox for allele-specific (AS) epigenomics analysis, which incorporates allelic variation data within existing
    resources, allowing for the identification of significant associations between epigenetic modifications and specific allelic variants
    in human and mouse cells. ALEA provides a customizable pipeline of command line tools for AS analysis of next-generation sequencing data
    (ChIP-seq, RNA-seq, etc.) that takes the raw sequencing data and produces separate allelic tracks ready to be viewed on genome browsers.
    ALEA takes advantage of the available genomic resources for human (The 1000 Genomes Project Consortium) and mouse (The Mouse Genome Project)
    to reconstruct diploid in silico genomes for human samples or hybrid mouse samples under study. Then, for each accompanying ChIP-seq or
    RNA-seq dataset, ALEA generates two wig files from short reads aligned differentially to each haplotype.
    This pipeline has been validated using human and hybrid mouse ChIP-seq and RNA-seq data (See Test Data section).
  doap:homepage: "http://www.bcgsc.ca/platform/bioinfo/software/alea"
  dcat:downloadURL: "ftp://ftp.bcgsc.ca/supplementary/ALEA/files/alea.1.2.2.tar.gz"
  doap:release:
  - class: doap:Version
    doap:revision: "1.2.2"
  doap:license: "AFL"
  doap:category: "commandline tool"
  doap:programming-language: "JAVA"
  foaf:publications:
  - id: urn:pmid:24371156
    foaf:title: >
      Hamid Younesy, Torsten Moller, Alireza Heravi-Moussavi, Jeffrey B. Cheng,
      Joseph F. Costello, Matthew C. Lorincz, Mohammad M. Karimi, Steven J. M. Jones
      ALEA: a toolbox for allele-specific epigenomics analysis Bioinformatics (2014) 30 (8): 1172-1174.
      doi: 10.1093/bioinformatics/btt744
    foaf:homepage: "http://bioinformatics.oxfordjournals.org/content/30/8/1172.long"
  doap:developer:
  - class: foaf:Organization
    foaf:name: "Canada's Michael Smith Genome Sciences Centre, BC Cancer Agency, Vancouver, British Columbia, V5Z 4S6, Canada"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Mohammad Karimi"
      foaf:mbox: "mailto:mkarimi@bcgsc.ca"
      foaf:homepage: "http://www.bcgsc.ca/author/mkarimi"

description: |
  alea-phaseVCF.cwl is developed for CWL consortium

doap:name: "alea-phaseVCF.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/alea-phaseVCF.cwl"

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
  - $import: alea-docker.cwl
  - class: InlineJavascriptRequirement

inputs:
  - id: "#hapsDir"
    type: File
    description: |
      path to the directory containing the .haps files
    inputBinding:
      position: 2

  - id: "#unphased"
    type: File
    description: |
      path to the vcf file containing unphased SNPs and Indels
    inputBinding:
      position: 3

  - id: "#outputPrefix"
    type: string
    description: |
      output file prefix including the path but not the extension
    inputBinding:
      position: 3

outputs:
  - id: "#phasevcf"
    type: File
    description: "Creates the file outputPrefix.vcf.gz"
    outputBinding:
      glob: $(inputs.outputPrefix+".vcf.gz")

baseCommand: ["alea", "phaseVCF"]

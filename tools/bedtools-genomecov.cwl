#!/usr/bin/env cwl-runner

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  adms: http://www.w3.org/ns/adms#
  dcat: http://www.w3.org/ns/dcat#
  edam: http://edamontology.org/

$schemas:
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf
- http://edamontology.org/EDAM.owl

cwlVersion: v1.0
class: CommandLineTool

adms:includedAsset:
  doap:name: bedtools
  doap:description: |
    A software suite for the comparison, manipulation and annotation of genomic features in browser extensible data (BED) and general feature format (GFF) format.
    BEDTools also supports the comparison of sequence alignments in BAM format to both BED and GFF features.
    The tools are extremely efficient and allow the user to compare large datasets (e.g. next-generation sequencing data) with both public and custom genome annotation tracks.
    BEDTools can be combined with one another as well as with standard UNIX commands, thus facilitating routine genomics tasks as well as pipelines that can quickly answer intricate questions of large genomic datasets.
  doap:homepage: http://bedtools.readthedocs.org
  doap:repository:
  - class: doap:GitRepository
    doap:location: https://github.com/arq5x/bedtools2
  doap:release:
  - class: doap:Version
    doap:revision: 2.25.0
  doap:license: GPLv2
  doap:category: commandline tool
  doap:programming-language: C++
  foaf:publications:
  - id: urn:pmid:20110278
    foaf:title: 'Aaron R. Quinlan, Ira M. Hall (2010) BEDTools: a flexible suite of
      utilities for comparing genomic features. Bioinformatics, 26(6) 841-842, http://dx.doi.org/10.1093/bioinformatics/btq033'
    foaf:homepage: http://bioinformatics.oxfordjournals.org/content/26/6/841
  doap:maintainer:
  - class: foaf:Person
    foaf:name: Aaron R. Quinlan
    foaf:mbox: aaronquinlan at gmail.com
    dct:isPartOf:
    - class: foaf:Organization
      foaf:name: Department of Biochemistry and Molecular Genetics, University of
        Virginia School of Medicine
    - class: foaf:Organization
      foaf:name: Center for Public Health Genomics, University of Virginia, Charlottesville,
        VA 22908, USA
doap:name: bedtools-genomecov.cwl
dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/bedtools-genomecov.cwl
doap:maintainer:
- class: foaf:Organization
  foaf:name: Barski Lab, Cincinnati Children's Hospital Medical Center
  foaf:member:
  - class: foaf:Person
    id: http://orcid.org/0000-0001-9102-5681
    foaf:openid: http://orcid.org/0000-0001-9102-5681
    foaf:name: Andrey Kartashov
    foaf:mbox: mailto:Andrey.Kartashov@cchmc.org
requirements:
- $import: bedtools-docker.yml
- class: InlineJavascriptRequirement

inputs:
  genomecoverageout: string
  scale:
    type: float?
    inputBinding:
      position: 4
      prefix: -scale

    doc: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
  max:
    type: int?
    inputBinding:
      position: 4
      prefix: -max
    doc: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
  m5:
    type: boolean?
    inputBinding:
      position: 4
      prefix: '-5'
    doc: |
      Calculate coverage of 5" positions (instead of entire interval).
  dept:
    type:
      name: JustDepts
      type: enum
      symbols: [-bg, -bga, -d]
    inputBinding:
      position: 4

  pairchip:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -pc
    doc: pair-end chip seq experiment
  dz:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -dz
    doc: |
      Report the depth at each genome position (with zero-based coordinates).
      Reports only non-zero positions.
      Default behavior is to report a histogram.
  split:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -split
    doc: |
      reat "split" BAM or BED12 entries as distinct BED intervals.
      when computing coverage.
      For BAM files, this uses the CIGAR "N" and "D" operations
      to infer the blocks for computing coverage.
      For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
      fields (i.e., columns 10,11,12).
  m3:
    type: boolean?
    inputBinding:
      position: 4
      prefix: '-3'
    doc: |
      Calculate coverage of 3" positions (instead of entire interval).
  fragmentsize:
    type: int?
    inputBinding:
      position: 4
      prefix: -fs
    doc: fixed fragment size
  input:
    type: File
    format: edam:format_2572
    inputBinding:
      position: 1
      valueFrom: |
        ${
          var prefix = ((/.*\.bam$/i).test(inputs.input.path))?'-ibam':'-i';
          return [prefix,inputs.input.path];
        }
    secondaryFiles: |
      ${
       if ((/.*\.bam$/i).test(inputs.input.path))
          return {"path": inputs.input.path+".bai", "class": "File"};
       return [];
      }

    doc: |
      The input file can be in BAM format
          (Note: BAM _must_ be sorted by position)
      or <bed/gff/vcf>
  strand:
    type: string?
    inputBinding:
      position: 4
      prefix: -strand
    doc: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6).
      - (STRING): can be + or -
  genomeFile:
    type: File
    inputBinding:
      position: 2
      prefix: -g
    doc: Input genome file.
outputs:
  genomecoverage:
    type: File
    outputBinding:
      glob: $(inputs.genomecoverageout)

    doc: The file containing the genome coverage
stdout: $(inputs.genomecoverageout)

baseCommand: [bedtools, genomecov]
doc: |
  bedtools-genomecov.cwl is developed for CWL consortium

  Original tool usage:
      Tool:    bedtools genomecov (aka genomeCoverageBed)
      Sources: https://github.com/arq5x/bedtools2
      Summary: Compute the coverage of a feature file among a genome.
      Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>


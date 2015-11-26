#!/usr/bin/env cwl-runner

"@context":
  "foaf": "http://xmlns.com/foaf/0.1/"
  "doap": "http://usefulinc.com/ns/doap"
  "adms": "http://purl.org/adms/"
  "admssw": "http://purl.org/adms/sw/"
  
adms:Asset:
  admssw:SoftwareProject:
    doap:name: "bedtools"
    doap:description: >
      A software suite for the comparison, manipulation and annotation of genomic features in browser extensible data (BED) and general feature format (GFF) format.
      BEDTools also supports the comparison of sequence alignments in BAM format to both BED and GFF features.
      The tools are extremely efficient and allow the user to compare large datasets (e.g. next-generation sequencing data) with both public and custom genome annotation tracks.
      BEDTools can be combined with one another as well as with standard UNIX commands, thus facilitating routine genomics tasks as well as pipelines that can quickly answer intricate questions of large genomic datasets.
    doap:homepage: "http://bedtools.readthedocs.org"
    doap:release:
      - doap:revision: "2.25.0"
    doap:license: "GPLv2"
    doap:category: "commandline tool"
    doap:programming-language: "C++"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/arq5x/bedtools2"
    foaf:Organization:
      - foaf:name: "Department of Biochemistry and Molecular Genetics, University of Virginia School of Medicine"
      - foaf:name: "Center for Public Health Genomics, University of Virginia, Charlottesville, VA 22908, USA"
    foaf:publications:
      - foaf:title: "(Quinlan and Hall, 2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics."
        foaf:homepage: "http://www.ncbi.nlm.nih.gov/pubmed/20110278"
    doap:maintainer:
      - foaf:Person:
        foaf:name: "Aaron R. Quinlan"
        foaf:mbox: "aaronquinlan at gmail.com"
  adms:AssetDistribution:
    doap:name: "bedtools-genomecov.cwl"
    doap:description: "Developed for CWL consortium http://commonwl.org/"
    doap:specification: "http://common-workflow-language.github.io/draft-3/"
    doap:release: "cwl:draft-3.dev2"
    doap:homepage: "http://commonwl.org/"
    doap:location: "https://github.com/common-workflow-language/workflows/blob/master/tools/bedtools-genomecov.cwl"
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
  Tool:    bedtools genomecov (aka genomeCoverageBed)
  Sources: https://github.com/arq5x/bedtools2
  Summary: Compute the coverage of a feature file among a genome.
  Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>

requirements:
  - "@import": envvar-global.cwl
  - "@import": bedtools-docker.cwl
  - class: InlineJavascriptRequirement


inputs:
  - id: "#input"
    type: File
    description: |
      The input file can be in BAM format
          (Note: BAM _must_ be sorted by position)
      or <bed/gff/vcf>
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

  - id: "#genomeFile"
    type: File
    description:
      Input genome file.
    inputBinding:
      position: 2
      prefix: "-g"

  - id: "#dept"
    type:
      name: "JustDepts"
      type: enum
      symbols: ["-bg","-bga","-d"]
    inputBinding:
      position: 4

  - id: "#scale"
    type: ["null",float ]
    description: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
    inputBinding:
      position: 4
      prefix: -scale

  - id: "#dz"
    type: ["null",boolean]
    description: |
      Report the depth at each genome position (with zero-based coordinates).
      Reports only non-zero positions.
      Default behavior is to report a histogram.
    inputBinding:
      position: 4
      prefix: "-dz"

  - id: "#split"
    type: ["null",boolean]
    description: |
      reat "split" BAM or BED12 entries as distinct BED intervals.
      when computing coverage.
      For BAM files, this uses the CIGAR "N" and "D" operations
      to infer the blocks for computing coverage.
      For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
      fields (i.e., columns 10,11,12).
    inputBinding:
      position: 4
      prefix: "-split"

  - id: "#strand"
    type: ["null", string]
    description: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6).
      - (STRING): can be + or -
    inputBinding:
      position: 4
      prefix: "-strand"

  - id: "#pairchip"
    type: ["null", boolean]
    description: "pair-end chip seq experiment"
    inputBinding:
      position: 4
      prefix: "-pc"

  - id: "#fragmentsize"
    type: ["null", int]
    description: "fixed fragment size"
    inputBinding:
      position: 4
      prefix: "-fs"

  - id: "#max"
    type: ["null",int]
    description: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
    inputBinding:
      position: 4
      prefix: "-max"

  - id: "#m5"
    type: ["null",boolean]
    description: |
      Calculate coverage of 5" positions (instead of entire interval).
    inputBinding:
      position: 4
      prefix: "-5"

  - id: "#m3"
    type: ["null",boolean]
    description: |
      Calculate coverage of 3" positions (instead of entire interval).
    inputBinding:
      position: 4
      prefix: "-3"

  - id: "#genomecoverageout"
    type: string

outputs:
  - id: "#genomecoverage"
    type: File
    description: "The file containing the genome coverage"
    outputBinding:
      glob: $(inputs.genomecoverageout)

stdout: $(inputs.genomecoverageout)

baseCommand: ["bedtools", "genomecov"]

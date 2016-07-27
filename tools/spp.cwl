#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool

requirements:
  - $import: spp-docker.yml
  - class: InlineJavascriptRequirement

inputs:
# ---------------------------------- #
# ------  MANDATORY ARGUMENTS  ----- #
# ---------------------------------- #
  - id: input_bam
    type: File
    description: "<ChIP_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped)(FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)"
    inputBinding:
      prefix: "-c="
      separate: false
# ------------------------------------------ #
# -- MANDATORY ARGUMENTS FOR PEAK CALLING -- #
# ------------------------------------------ #
  - id: control_bam
    type:
      - 'null'
      - File
    description: "<Input_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)"
    inputBinding:
      prefix: "-i="
      separate: false
# ---------------------------------- #
# -------  OPTIONAL ARGUMENTS  ----- #
# ---------------------------------- #
  - id: nthreads
    type:
      - 'null'
      - int
    description: "-p=<nodes> , number of parallel processing nodes, default=0"
    inputBinding:
      prefix: "-p="
      separate: false
  - id: s
    type:
      - 'null'
      - string
    description: "-s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated, default=-500:5:1500"
    inputBinding:
      prefix: "-s="
      separate: false
  - id: speak
    type:
      - 'null'
      - string
    description: "-speak=<strPeak>, user-defined cross-correlation peak strandshift"
    inputBinding:
      prefix: "-speak="
      separate: false
  - id: x
    type:
      - 'null'
      - string
    description: "-x=<min>:<max>, strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10)"
    inputBinding:
      prefix: "-x="
      separate: false
  - id: fdr
    type:
      - 'null'
      - float
    description: "-fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling"
    inputBinding:
      prefix: "-fdr="
      separate: false
  - id: npeak
    type:
      - 'null'
      - int
    description: "-npeak=<numPeaks>, threshold on number of peaks to call"
    inputBinding:
      prefix: "-npeak="
      separate: false
  - id: filtchr
    type:
      - 'null'
      - string
    description: "-filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name"
    inputBinding:
      prefix: "-filtchr="
      separate: false
  - id: savp
    type:
      - 'null'
      - boolean
    description: "save cross-correlation plot"
    inputBinding:
      valueFrom: ${ if (self) return "-savp=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.pdf"; return null}
  - id: savn
    type:
      - 'null'
      - boolean
    description: "-savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed width peaks)"
    inputBinding:
      valueFrom: ${ if (self) return "-savn=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.narrowPeak"; return null}
  - id: savr
    type:
      - 'null'
      - boolean
    description: "-savr=<regionpeakfilename> OR -savr RegionPeak file name (variable width peaks with regions of enrichment)"
    inputBinding:
      valueFrom: ${ if (self) return "-savr=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.regionPeak"; return null}
  - id: savd
    type:
      - 'null'
      - boolean
    description: "-savd=<rdatafile> OR -savd, save Rdata file"
    inputBinding:
      valueFrom: ${ if (self) return "-savd=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.Rdata"; return null}
  - id: rf
    type: boolean
    default: true
    description: "overwrite (force remove) output files in case they exist. Default: true"
    inputBinding:
      prefix: "-rf"

outputs:
  - id: output_spp_cross_corr
    type: File
    description: "peakshift/phantomPeak results summary file"
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.txt")
  - id: output_spp_cross_corr_plot
    type: ['null', File]
    description: "peakshift/phantomPeak results summary plot"
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.pdf")
  - id: output_spp_narrow_peak
    type: ['null', File]
    description: "narrowPeak output file"
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.narrowPeak")
  - id: output_spp_region_peak
    type: ['null', File]
    description: "regionPeak output file"
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.regionPeak")
  - id: output_spp_rdata
    type: ['null', File]
    description: "Rdata file from the run_spp.R run"
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.Rdata")

baseCommand: run_spp.R
arguments:
  - valueFrom: $("-out=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.txt")

#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: spp-docker.yml
- class: InlineJavascriptRequirement

inputs:
# ---------------------------------- #
# ------  MANDATORY ARGUMENTS  ----- #
# ---------------------------------- #
  nthreads:
# ---------------------------------- #
# -------  OPTIONAL ARGUMENTS  ----- #
# ---------------------------------- #
    type: int?
    inputBinding:
      prefix: -p=
      separate: false
    doc: -p=<nodes> , number of parallel processing nodes, default=0
  filtchr:
    type: string?
    inputBinding:
      prefix: -filtchr=
      separate: false
    doc: -filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific
      chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their
      name
  savr:
    type: boolean?
    inputBinding:
      valueFrom: ${ if (self) return "-savr=" + inputs.input_bam.path.replace(/^.*[\\\/]/,
        "").replace(/\.[^/.]+$/, "") + ".spp.regionPeak"; return null}
    doc: -savr=<regionpeakfilename> OR -savr RegionPeak file name (variable width
      peaks with regions of enrichment)
  rf:
    type: boolean
    default: true
    inputBinding:
      prefix: -rf
    doc: 'overwrite (force remove) output files in case they exist. Default: true'
  npeak:
    type: int?
    inputBinding:
      prefix: -npeak=
      separate: false
    doc: -npeak=<numPeaks>, threshold on number of peaks to call
  savn:
    type: boolean?
    inputBinding:
      valueFrom: ${ if (self) return "-savn=" + inputs.input_bam.path.replace(/^.*[\\\/]/,
        "").replace(/\.[^/.]+$/, "") + ".spp.narrowPeak"; return null}
    doc: -savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed width peaks)
  s:
    type: string?
    inputBinding:
      prefix: -s=
      separate: false
    doc: -s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated,
      default=-500:5:1500
  fdr:
    type: float?
    inputBinding:
      prefix: -fdr=
      separate: false
    doc: -fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling
  savp:
    type: boolean?
    inputBinding:
      valueFrom: ${ if (self) return "-savp=" + inputs.input_bam.path.replace(/^.*[\\\/]/,
        "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.pdf"; return null}
    doc: save cross-correlation plot
  x:
    type: string?
    inputBinding:
      prefix: -x=
      separate: false
    doc: -x=<min>:<max>, strand shifts to exclude (This is mainly to avoid region
      around phantom peak) default=10:(readlen+10)
  control_bam:
# ------------------------------------------ #
# -- MANDATORY ARGUMENTS FOR PEAK CALLING -- #
# ------------------------------------------ #
    type: File?
    inputBinding:
      prefix: -i=
      separate: false
    doc: <Input_alignFile>, full path and name (or URL) of tagAlign/BAM file (can
      be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
  savd:
    type: boolean?
    inputBinding:
      valueFrom: ${ if (self) return "-savd=" + inputs.input_bam.path.replace(/^.*[\\\/]/,
        "").replace(/\.[^/.]+$/, "") + ".spp.Rdata"; return null}
    doc: -savd=<rdatafile> OR -savd, save Rdata file
  input_bam:
    type: File
    inputBinding:
      prefix: -c=
      separate: false
    doc: <ChIP_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be
      gzipped)(FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
  speak:
    type: string?
    inputBinding:
      prefix: -speak=
      separate: false
    doc: -speak=<strPeak>, user-defined cross-correlation peak strandshift
outputs:
  output_spp_cross_corr:
    type: File
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/,
        "") + ".spp_cross_corr.txt")
    doc: peakshift/phantomPeak results summary file
  output_spp_rdata:
    type: File?
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/,
        "") + ".spp.Rdata")

    doc: Rdata file from the run_spp.R run
  output_spp_narrow_peak:
    type: File?
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/,
        "") + ".spp.narrowPeak")
    doc: narrowPeak output file
  output_spp_region_peak:
    type: File?
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/,
        "") + ".spp.regionPeak")
    doc: regionPeak output file
  output_spp_cross_corr_plot:
    type: File?
    outputBinding:
      glob: $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/,
        "") + ".spp_cross_corr.pdf")
    doc: peakshift/phantomPeak results summary plot
baseCommand: run_spp.R
arguments:
- valueFrom: $("-out=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/,
    "") + ".spp_cross_corr.txt")


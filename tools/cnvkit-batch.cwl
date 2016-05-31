#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'batch']

requirements:
  - class: InlineJavascriptRequirement

description: |
  Run the complete CNVkit pipeline on one or more BAM files.

inputs:


- id: bam_files
  type:
  - "null"
  - type: array
    items: string

  description: Mapped sequence reads (.bam)
  inputBinding:
    position: 1


- id: male_reference
  type: ["null", boolean]
  default: null
  description: Use or assume a male reference (i.e. female samples will have +1
                log-CNR of chrX; otherwise male samples would have -1 chrX).
  inputBinding:
    position: 2
    prefix: --male-reference 

- id: count_reads
  type: ["null", boolean]
  default: null
  description: Get read depths by counting read midpoints within each bin.
                (An alternative algorithm).
  inputBinding:
    position: 3
    prefix: --count-reads 

- id: processes
  type: ["null", int]
  default: 1
  description: Number of subprocesses used to running each of the BAM files in
                parallel. Give 0 or a negative value to use the maximum number
                of available CPUs. [Default - process each BAM in serial]
  inputBinding:
    position: 4
    prefix: --processes 

- id: rlibpath
  type: ["null", string]
  description: Path to an alternative site-library to use for R packages.
  inputBinding:
    position: 5
    prefix: --rlibpath 

- id: normal
  type:
  - "null"
  - type: array
    items: string

  description: Normal samples (.bam) to construct the pooled reference.
                If this option is used but no files are given, a "flat"
                reference will be built.
  inputBinding:
    position: 6
    prefix: --normal 

- id: fasta
  type: ["null", string]
  description: Reference genome, FASTA format (e.g. UCSC hg19.fa)
  inputBinding:
    position: 7
    prefix: --fasta 

- id: targets
  type: ["null", string]
  description: Target intervals (.bed or .list)
  inputBinding:
    position: 8
    prefix: --targets 

- id: antitargets
  type: ["null", string]
  description: Antitarget intervals (.bed or .list)
  inputBinding:
    position: 9
    prefix: --antitargets 

- id: annotate
  type: ["null", string]
  description: UCSC refFlat.txt or ensFlat.txt file for the reference genome.
                Pull gene names from this file and assign them to the target
                regions.
  inputBinding:
    position: 10
    prefix: --annotate 

- id: short_names
  type: ["null", boolean]
  default: null
  description: Reduce multi-accession bait labels to be short and consistent.
  inputBinding:
    position: 11
    prefix: --short-names 

- id: split
  type: ["null", boolean]
  default: null
  description: Split large tiled intervals into smaller, consecutive targets.
  inputBinding:
    position: 12
    prefix: --split 

- id: target_avg_size
  type: ["null", int]
  description: Average size of split target bins (results are approximate).
  inputBinding:
    position: 13
    prefix: --target-avg-size 

- id: access
  type: ["null", string]
  description: Regions of accessible sequence on chromosomes (.bed), as
                output by the 'access' command.
  inputBinding:
    position: 14
    prefix: --access 

- id: antitarget_avg_size
  type: ["null", int]
  description: Average size of antitarget bins (results are approximate).
  inputBinding:
    position: 15
    prefix: --antitarget-avg-size 

- id: antitarget_min_size
  type: ["null", int]
  description: Minimum size of antitarget bins (smaller regions are dropped).
  inputBinding:
    position: 16
    prefix: --antitarget-min-size 

- id: output_reference
  type: ["null", string]
  description: Output filename/path for the new reference file being created.
                (If given, ignores the -o/--output-dir option and will write the
                file to the given path. Otherwise, "reference.cnn" will be
                created in the current directory or specified output directory.)
                
  inputBinding:
    position: 17
    prefix: --output-reference 

- id: reference
  type: ["null", string]
  description: Copy number reference file (.cnn).
  inputBinding:
    position: 18
    prefix: --reference 

- id: output_dir
  type: ["null", string]
  default: .
  description: Output directory.
  inputBinding:
    position: 19
    prefix: --output-dir 

- id: scatter
  type: ["null", boolean]
  default: null
  description: Create a whole-genome copy ratio profile as a PDF scatter plot.
  inputBinding:
    position: 20
    prefix: --scatter 

- id: diagram
  type: ["null", boolean]
  default: null
  description: Create a diagram of copy ratios on chromosomes as a PDF.
  inputBinding:
    position: 21
    prefix: --diagram 
outputs:
    []

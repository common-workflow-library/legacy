#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'scatter']

requirements:
  - class: InlineJavascriptRequirement

description: |
  Plot probe log2 coverages and segmentation calls together.

inputs:
  

- id: filename
  type: ["null", string]
  description: Processed bin-level copy ratios (*.cnr), the output
                of the 'fix' sub-command.
  inputBinding:
    position: 1

- id: segment
  type: ["null", string]
  description: Segmentation calls (.cns), the output of the 'segment' command.
  inputBinding:
    prefix: --segment 

- id: chromosome
  type: ["null", string]
  description: Chromosome (e.g. 'chr1') or chromosomal range (e.g.
                'chr1 -2333000-2444000') to display. If a range is given,
                all targeted genes in this range will be shown, unless
                '--gene'/'-g' is already given.
  inputBinding:
    prefix: --chromosome 

- id: gene
  type: ["null", string]
  description: Name of gene or genes (comma-separated) to display.
  inputBinding:
    prefix: --gene 

- id: range_list
  type: ["null", string]
  description: File listing the chromosomal ranges to display, as BED, interval
                list or "chr -start-end" text. Creates focal plots similar to
                -c/--chromosome for each listed region, combined into a
                multi-page PDF.  The output filename must also be
                specified (-o/--output).
  inputBinding:
    prefix: --range-list 

- id: sample_id
  type: ["null", string]
  description: Specify the name of the sample in the VCF to use for b-allele
                frequency extraction and to show in plot title.
  inputBinding:
    prefix: --sample-id 

- id: normal_id
  type: ["null", string]
  description: Corresponding normal sample ID in the input VCF.
  inputBinding:
    prefix: --normal-id 

- id: background_marker
  type: ["null", string]
  description: Plot antitargets with this symbol, in zoomed/selected regions.
                [Default - same as targets]
  inputBinding:
    prefix: --background-marker 

- id: trend
  type: ["null", boolean]
  default: null
  description: Draw a smoothed local trendline on the scatter plot.
  inputBinding:
    prefix: --trend 

- id: vcf
  type: ["null", string]
  description: VCF file name containing variants to plot for SNV allele
                frequencies.
  inputBinding:
    prefix: --vcf 

- id: min_variant_depth
  type: ["null", int]
  default: 20
  description: Minimum read depth for a SNV to be displayed in the b-allele
                frequency plot. [Default - %(default)s]
  inputBinding:
    prefix: --min-variant-depth 

- id: width
  type: ["null", float]
  default: 1000000.0
  description: Width of margin to show around the selected gene or region
                on the chromosome (use with --gene or --region).
                [Default - %(default)d]
  inputBinding:
    prefix: --width 

- id: y_min
  type: ["null", float]
  description: y-axis lower limit.
  inputBinding:
    prefix: --y-min 

- id: y_max
  type: ["null", float]
  description: y-axis upper limit.
  inputBinding:
    prefix: --y-max 

- id: output
  type: ["null", string]
  description: Output table file name.
  inputBinding:
    prefix: --output 

outputs:
    []

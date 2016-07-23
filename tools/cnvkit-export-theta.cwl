#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2cwl ver. 0.2.8
# To generate again: $ cnvkit.py --generate_cwl_tool
# Help: $ cnvkit.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['cnvkit.py', 'export', 'theta']

description: |
  Convert segments to THetA2 input file format (*.input).

inputs:
  
  tumor_segment:
    type: string
  
    description: Tumor-sample segmentation file from CNVkit (.cns).
    inputBinding:
      position: 1

  normal_reference:
    type: ["null", string]
    description: Reference copy number profile (.cnn), or normal-sample bin-level log2 copy ratios (.cnr). [DEPRECATED]
    inputBinding:
      position: 2

  reference:
    type: ["null", string]
    description: Reference copy number profile (.cnn), or normal-sample bin-level log2 copy ratios (.cnr). Use if the tumor_segment input file does not contain a "weight" column.
    inputBinding:
      prefix: --reference 

  vcf:
    type: ["null", string]
    description: VCF file containing SNVs observed in both the tumor and normal samples. Tumor sample ID should match the `tumor_segment` filename or be specified with -i/--sample-id.
    inputBinding:
      prefix: --vcf 

  sample_id:
    type: ["null", string]
    description: Specify the name of the tumor sample in the VCF (given with -v/--vcf). [Default - taken the tumor_segment file name]
    inputBinding:
      prefix: --sample-id 

  normal_id:
    type: ["null", string]
    description: Corresponding normal sample ID in the input VCF.
    inputBinding:
      prefix: --normal-id 

  min_depth:
    type: ["null", int]
    default: 20
    description: Minimum read depth for a SNP in the VCF to be counted. [Default - %(default)s]
    inputBinding:
      prefix: --min-depth 

  output:
    type: ["null", string]
    description: Output file name.
    inputBinding:
      prefix: --output 


outputs:
    []

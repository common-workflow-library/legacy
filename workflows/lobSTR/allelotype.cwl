#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

inputs:
  noise_model:

    type: File
    inputBinding:
      prefix: --noise_model
      valueFrom: |
        ${ return {"path": self.path.match(/(.*)\.stepmodel/)[1], "class": "File"}; }
    secondaryFiles:
    - ^.stuttermodel

    doc: |
      File to read noise model parameters from (.stepmodel)
  bam:
    type: File
    inputBinding:
      prefix: --bam
    secondaryFiles:
    - .bai

    doc: |
      BAM file to analyze. Should have a unique read group and be sorted and indexed.
  reference:

    type: File
    inputBinding:
      prefix: --index-prefix
      valueFrom: |
        ${ return {"path": self.path.match(/(.*)ref\.fasta/)[1], "class": "File"}; }

    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .pac
    - .rbwt
    - .rpac
    - .rsa
    - .sa
    - ${return self.basename.replace(/(.*)ref\.fasta/, "$1chromsizes.tab");}
    - ${return self.basename.replace(/(.*)ref\.fasta/, "$1mergedref.bed");}
    - ${return self.basename.replace(/(.*)ref\.fasta/, "$1ref_map.tab");}

    doc: lobSTR's bwa reference files
  strinfo:
    type: File
    inputBinding:
      prefix: --strinfo
    doc: |
      File containing statistics for each STR.
  output_prefix:
    type: string
    inputBinding:
      prefix: --out
    doc: Prefix for output files. will output prefix.vcf and prefix.genotypes.tab
outputs:
  vcf_stats:
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.allelotype.stats')

  vcf:
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.vcf')
baseCommand: [allelotype, --command, classify]

arguments:
- --noweb
doc: Run lobSTR allelotype classifier.


#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: superenhancer.yml
- class: InlineJavascriptRequirement

inputs:
  inputbai:
    type: File
    doc: Input indexed BAI file
  stitch_distance:
    type: int?
    inputBinding:
      prefix: -s
    doc: 'max distance between two regios to stitch together. (Default: 12500)'
  controlbam:
    type: File?
    inputBinding:
      prefix: -c
    doc: .bam file used as control. Subtracted from the density of the ranking_Bam
  inputbam:
    type: File
    inputBinding:
      prefix: -r
    doc: Input sorted BAM file
  inputgff:
    type: File
    inputBinding:
      prefix: -i
    doc: Input GFF file
  genome:
    type: File
    inputBinding:
      prefix: -g
    doc: 'REFERENCE GENOME BUILD UCSC FILE: one of hg18, hg19, mm8, mm9, mm10 file
      for example hg18_refseq.ucsc'
  output_name:
    type: string
    inputBinding:
      prefix: -o
    doc: output folder name
  tss:
    type: int?
    inputBinding:
      prefix: -t
    doc: 'TSS_Exclusion: excludes regios contained within +/- this distance from TSS
      in order to account for promoter bias (Default: 0, recommended: 2500)'
  bams:
    type: File?
    inputBinding:
      prefix: -b
    doc: A comma separated list of additional bam files to map to
outputs:
  All_enhancers:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.output_name)


baseCommand: [python, /usr/local/bin/ROSE_main.py]
$schemas:
- http://schema.org/docs/schema_org_rdfa.html

$namespaces:
  s: http://schema.org/

s:mainEntity:
  class: s:SoftwareSourceCode
  s:name: superenhnacer
  s:about: ''
  PURPOSE: To create stitched enhancers, and to separate super-enhancers from typical
    enhancers using sequencing data (.bam) given a file of previously identified constituent
    enhancers (.gff)It makes use of the superenhancer script developed by Young Lab

  s:url: http://younglab.wi.mit.edu/super_enhancer_code.html

  s:codeRespository: https://github.com/bharath-cchmc/edited-Super-Enhancer

  s:license:
  - http://younglab.wi.mit.edu/ROSE/LICENSE.txt
  - https://opensource.org/licenses/MIT
  - https://opensource.org/licenses/BSD-3-Clause

  s:targetProduct:
    class: s:SoftwareApplication
    s:applicationCategory: CommandLine Tool
  s:programmingLanguage: Python
  s:publication:
  - class: s:ScholarlyArticle
    id: https://doi.org/10.1016/j.cell.2013.03.035
    s:name: Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles
      Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young Cell
    s:url: http://www.cell.com/abstract/S0092-8674(13)00392-9

  - class: s:ScholarlyArticle
    id: https://doi.org/10.1016/j.cell.2013.03.036
    s:name: Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando,
      Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young Cell
    s:url: http://www.cell.com/abstract/S0092-8674(13)00393-0

s:downloadUrl: https://bitbucket.org/bharath-cchmc/se-docker-cwl-trial
s:codeRepository: https://github.com/common-workflow-language/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:author:
  class: s:Person
  s:name: Bharath Manicka Vasagam
  s:email: mailto:Bharath.Manickavasagam@cchmc.org
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab

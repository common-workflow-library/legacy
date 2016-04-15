#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: samtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#]'
    type: boolean
    description: ']'
    inputBinding:
      position: 4
  - id: '#in1bam'
    type: File
    description: in1.bam
    inputBinding:
      position: 2
  - id: '#6'
    type:
      - 'null'
      - boolean
    description: |
      --illumina1.3+      quality is in the Illumina-1.3+ encoding
    inputBinding:
      position: 1
      prefix: '-6'
  - id: '#A'
    type:
      - 'null'
      - boolean
    description: |
      --count-orphans     do not discard anomalous read pairs
    inputBinding:
      position: 1
      prefix: '-A'
  - id: '#b'
    type:
      - 'null'
      - File
    description: |
      --bam-list FILE     list of input BAM filenames, one per line
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#B'
    type:
      - 'null'
      - boolean
    description: |
      --no-BAQ            disable BAQ (per-Base Alignment Quality)
    inputBinding:
      position: 1
      prefix: '-B'
  - id: '#C'
    type:
      - 'null'
      - int
    description: |
      --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
    inputBinding:
      position: 1
      prefix: '-C'
  - id: '#d'
    type:
      - 'null'
      - int
    description: >
      --max-depth INT     max per-BAM depth; avoids excessive memory usage [250]
    inputBinding:
      position: 1
      prefix: '-d'
  - id: '#E'
    type:
      - 'null'
      - boolean
    description: |
      --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs
    inputBinding:
      position: 1
      prefix: '-E'
  - id: '#f'
    type:
      - 'null'
      - File
    description: |
      --fasta-ref FILE    faidx indexed reference sequence file
    inputBinding:
      position: 1
      prefix: '-f'
  - id: '#G'
    type:
      - 'null'
      - File
    description: |
      --exclude-RG FILE   exclude read groups listed in FILE
    inputBinding:
      position: 1
      prefix: '-G'
  - id: '#l'
    type:
      - 'null'
      - boolean
    description: |
      --positions FILE    skip unlisted positions (chr pos) or regions (BED)
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#q'
    type:
      - 'null'
      - int
    description: |
      --min-MQ INT        skip alignments with mapQ smaller than INT [0]
    inputBinding:
      position: 1
      prefix: '-q'
  - id: '#Q'
    type:
      - 'null'
      - int
    description: |
      --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
    inputBinding:
      position: 1
      prefix: '-Q'
  - id: '#r'
    type:
      - 'null'
      - boolean
    description: |
      --region REG        region in which pileup is generated
    inputBinding:
      position: 1
      prefix: '-r'
  - id: '#R'
    type:
      - 'null'
      - int
    description: |
      --ignore-RG         ignore RG tags (one BAM = one sample)
      --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset []
      --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set
      [UNMAP,SECONDARY,QCFAIL,DUP]
    inputBinding:
      position: 1
      prefix: '-R'
  - id: '#x'
    type:
      - 'null'
      - boolean
    description: |
      --ignore-overlaps   disable read-pair overlap detection
      Output options:
    inputBinding:
      position: 1
      prefix: '-x'
  - id: '#o'
    type:
      - 'null'
      - File
    description: |
      --output FILE       write output to FILE [standard output]
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#g'
    type:
      - 'null'
      - boolean
    description: |
      --BCF               generate genotype likelihoods in BCF format
    inputBinding:
      position: 1
      prefix: '-g'
  - id: '#v'
    type:
      - 'null'
      - boolean
    description: |
      --VCF               generate genotype likelihoods in VCF format
      Output options for mpileup format (without -g/-v):
    inputBinding:
      position: 1
      prefix: '-v'
  - id: '#O'
    type:
      - 'null'
      - boolean
    description: |
      --output-BP         output base positions on reads
    inputBinding:
      position: 1
      prefix: '-O'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      --output-MQ         output mapping quality
      Output options for genotype likelihoods (when -g/-v is used):
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#t'
    type:
      - 'null'
      - boolean
    description: |
      --output-tags LIST  optional tags to output: DP,DPR,DV,DP4,INFO/DPR,SP []
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#u'
    type:
      - 'null'
      - boolean
    description: |
      --uncompressed      generate uncompressed VCF/BCF output
      SNP/INDEL genotype likelihoods options (effective with -g/-v):
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#e'
    type:
      - 'null'
      - int
    description: |
      --ext-prob INT      Phred-scaled gap extension seq error probability [20]
    inputBinding:
      position: 1
      prefix: '-e'
  - id: '#F'
    type:
      - 'null'
      - float
    description: |
      --gap-frac FLOAT    minimum fraction of gapped reads [0.002]
    inputBinding:
      position: 1
      prefix: '-F'
  - id: '#h'
    type:
      - 'null'
      - int
    description: |
      --tandem-qual INT   coefficient for homopolymer errors [100]
    inputBinding:
      position: 1
      prefix: '-h'
  - id: '#I'
    type:
      - 'null'
      - boolean
    description: |
      --skip-indels       do not perform indel calling
    inputBinding:
      position: 1
      prefix: '-I'
  - id: '#L'
    type:
      - 'null'
      - int
    description: |
      --max-idepth INT    maximum per-sample depth for INDEL calling [250]
    inputBinding:
      position: 1
      prefix: '-L'
  - id: '#m'
    type:
      - 'null'
      - int
    description: |
      --min-ireads INT    minimum number gapped reads for indel candidates [1]
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#o'
    type:
      - 'null'
      - int
    description: |
      --open-prob INT     Phred-scaled gap open seq error probability [40]
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#p'
    type:
      - 'null'
      - boolean
    description: |
      --per-sample-mF     apply -m and -F per-sample for increased sensitivity
    inputBinding:
      position: 1
      prefix: '-p'
  - id: '#P'
    type:
      - 'null'
      - string
    description: |
      --platforms STR     comma separated list of platforms for indels [all]
      Notes: Assuming diploid individuals.
    inputBinding:
      position: 1
      prefix: '-P'
outputs:
  - id: '#stdoutfile'
    type: File
    outputBinding:
      glob:
        engine: 'cwl:JsonPointer'
        script: /job/stdoutfile
stdout:
  engine: 'cwl:JsonPointer'
  script: /job/stdoutfile
baseCommand:
  - samtools
  - mpileup
description: |-
  Usage: samtools mpileup [options] in1.bam [in2.bam [...]]

  Input options:
    -6, --illumina1.3+      quality is in the Illumina-1.3+ encoding
    -A, --count-orphans     do not discard anomalous read pairs
    -b, --bam-list FILE     list of input BAM filenames, one per line
    -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)
    -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
    -d, --max-depth INT     max per-BAM depth; avoids excessive memory usage [250]
    -E, --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs
    -f, --fasta-ref FILE    faidx indexed reference sequence file
    -G, --exclude-RG FILE   exclude read groups listed in FILE
    -l, --positions FILE    skip unlisted positions (chr pos) or regions (BED)
    -q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
    -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
    -r, --region REG        region in which pileup is generated
    -R, --ignore-RG         ignore RG tags (one BAM = one sample)
    --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset []
    --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set
                                              [UNMAP,SECONDARY,QCFAIL,DUP]
    -x, --ignore-overlaps   disable read-pair overlap detection

  Output options:
    -o, --output FILE       write output to FILE [standard output]
    -g, --BCF               generate genotype likelihoods in BCF format
    -v, --VCF               generate genotype likelihoods in VCF format

  Output options for mpileup format (without -g/-v):
    -O, --output-BP         output base positions on reads
    -s, --output-MQ         output mapping quality

  Output options for genotype likelihoods (when -g/-v is used):
    -t, --output-tags LIST  optional tags to output: DP,DPR,DV,DP4,INFO/DPR,SP []
    -u, --uncompressed      generate uncompressed VCF/BCF output

  SNP/INDEL genotype likelihoods options (effective with -g/-v):
    -e, --ext-prob INT      Phred-scaled gap extension seq error probability [20]
    -F, --gap-frac FLOAT    minimum fraction of gapped reads [0.002]
    -h, --tandem-qual INT   coefficient for homopolymer errors [100]
    -I, --skip-indels       do not perform indel calling
    -L, --max-idepth INT    maximum per-sample depth for INDEL calling [250]
    -m, --min-ireads INT    minimum number gapped reads for indel candidates [1]
    -o, --open-prob INT     Phred-scaled gap open seq error probability [40]
    -p, --per-sample-mF     apply -m and -F per-sample for increased sensitivity
    -P, --platforms STR     comma separated list of platforms for indels [all]

  Notes: Assuming diploid individuals.


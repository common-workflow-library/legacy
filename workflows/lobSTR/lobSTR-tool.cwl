#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: lobSTR
description: "lobSTR is a tool for profiling Short Tandem Repeats (STRs) from high throughput sequencing data."
#"doap:homepage": http://melissagymrek.com/lobstr-code/index.html
#"doap:license": http://spdx.org/licenses/GPL-3.0
requirements:
  - class: InlineJavascriptRequirement

inputs:

  - id: files
    type:
      - "null"
      - type: array
        items: File
    description: List of files of single ends in fasta, fastq, or BAM files to align. Also use this parameter to specify paired-end BAM files to align.
    #format: ["edam:format_1929", "edam:format_1930", "edam:format_2572"]
    inputBinding:
      prefix: "--files"
      position: 2
      itemSeparator: ','

  - id: p1
    type:
      - "null"
      - type: array
        items: File
    description: list of files containing the first end of paired end reads in fasta or fastq format
    #format: ["edam:format_1929", "edam:format_1930"]
    inputBinding:
      prefix: "--p1"
      position: 2
      itemSeparator: ','

  - id: p2
    type:
      - "null"
      - type: array
        items: File
    description: list of files containing the second end of paired end reads in fasta or fastq format
    #format: ["edam:format_1929", "edam:format_1930"]
    inputBinding:
      prefix: "--p2"
      position: 3
      itemSeparator: ','

  - id: output_prefix
    type: string
    description: prefix for output files. will output prefix.aligned.bam and prefix.aligned.stats
    inputBinding:
      prefix: "--out"

  - id: reference
    type: File
    description: "lobSTR's bwa reference files"
    #format: "edam:format_1929"
    inputBinding:
      prefix: "--index-prefix"
      valueFrom: |
          ${ return {"path": self.path.match(/(.*)ref\.fasta/)[1], "class": "File"}; }

    secondaryFiles:
      - ".amb"
      - ".ann"
      - ".bwt"
      - ".pac"
      - ".rbwt"
      - ".rpac"
      - ".rsa"
      - ".sa"
      - ${return self.path.replace(/(.*)ref\.fasta/, "$1chromsizes.tab");}
      - ${return self.path.replace(/(.*)ref\.fasta/, "$1mergedref.bed");}
      - ${return self.path.replace(/(.*)ref\.fasta/, "$1ref_map.tab");}

  - id: rg-sample
    type: string
    inputBinding:
      prefix: "--rg-sample"
    description: Use this in the read group SM tag

  - id: rg-lib
    type: string
    inputBinding:
      prefix: "--rg-lib"
    description: Use this in the read group LB tag

  - id: bampair
    type: ["null", boolean]
    inputBinding:
      prefix: "--bampair"
    description: |
      Reads are in BAM format and are paired-end. NOTE: BAM file MUST be sorted or collated by read name.

  - id: bwaq
    type: ["null", int]
    inputBinding:
      prefix: "--bwaq"
    description: trim reads based on quality score. Same as the -q parameter in BWA

  - id: oldillumina
    type: ["null", boolean]
    inputBinding:
      prefix: "--oldillumina"
    description: |
      specifies that base pair quality scores are reported in the old Phred
      format (Illumina 1.3+, Illumina 1.5+) where quality scores are given as
      Phred + 64 rather than Phred + 33

  - id: multi
    type: ["null", boolean]
    inputBinding:
      prefix: "--multi"
    description: Report reads mapping to multiple genomic locations

  - id: min-read-length
    type: int
    default: 45
    inputBinding:
      prefix: "--min-read-length"
    description: Don't process reads shorter than this

  - id: max-read-length
    type: int
    default: 1024
    inputBinding:
      prefix: "--min-read-length"
    description: Don't process reads longer than this

  - id: fft-window-size
    type: int
    default: 16
    inputBinding:
      prefix: "--fft-window-size"
    description: Size of sliding entropy window

  - id: fft-window-step
    type: int
    default: 4
    inputBinding:
      prefix: "--fft-window-step"
    description: Step size of the sliding window

  - id: entropy-threshold
    type: float
    default: 0.45
    inputBinding:
      prefix: "--entropy-threshold"
    description: Threshold score to call a window periodic

  - id: minflank
    type: int
    default: 8
    inputBinding:
      prefix: "--minflank"
    description: Minimum length of flanking region to try to align

  - id: maxflank
    type: int
    default: 100
    inputBinding:
      prefix: "--maxflank"
    description: Length to trim flanking regions to before aligning

  - id: max-diff-ref
    type: int
    default: 50
    inputBinding:
      prefix: "--max-diff-ref"
    description: Only report reads differing by at most this number of bp from the reference allele

  - id: extend
    type: int
    default: 1000
    inputBinding:
      prefix: "--extend"
    description: Number of bp the reference was extended when building the index. Must be same as --extend parameter used to run lobstr_index.py

  - id: maq
    type: int
    default: 100
    inputBinding:
      prefix: "--mapq"
    description: maximum allowed mapq score calculated as the sum of qualities at base mismatches

  - id: u
    type: ["null", boolean]
    inputBinding:
      prefix: "-u"
    description: only report reads different by an integer number of copy numbers from the reference allele

  - id: mismatch
    type: int
    default: 1
    inputBinding:
      prefix: "--mismatch"
    description: allowed edit distance in either flanking region

  - id: g
    type: int
    default: 1
    inputBinding:
      prefix: "-g"
    description: Maximum number of gap opens allowed in each flanking region

  - id: e
    type: int
    default: 1
    inputBinding:
      prefix: "-e"
    description: Maximum number of gap extends allowed in each flanking region

  - id: r
    type: float
    default: -1
    inputBinding:
      prefix: "-e"
    description: fraction of missing alignments given 2% uniform base error rate. Ignored if --mismatch is set

  - id: max-hits-quit-aln
    type: int
    default: 1000
    inputBinding:
      prefix: "--max-hits-quit-aln"
    description: Stop alignment search after this many hits found. Use -1 for no limit.

  - id: min-flank-allow-mismatch
    type: int
    default: 30
    inputBinding:
      prefix: "--min-flank-allow-mismatch"
    description: Don't allow mismatches if aligning flanking regions shorter than this

outputs:
  - id: bam
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.aligned.bam')

  - id: bam_stats
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.aligned.stats')

baseCommand: ["lobSTR"]
arguments:
  - "--verbose"
  - "--noweb"
  - valueFrom: |
        ${
          var r = /.*(\.fastq|\.fq)(\.gz)?$/i;
          var m;
          if (inputs['files']) {
            m = r.exec(inputs.files[0].path);
          } else {
            m = r.exec(inputs.p1[0].path);
          }
          if (m) {
            if (m[2]) {
              return ["--fastq", "--gzip"];
            } else {
              return "--fastq";
            }
          } else {
            return null;
          }
        }
  - valueFrom: |
        ${
          if (inputs['files'] && (/.*\.bam$/i).test(inputs['files'][0].path) && !inputs.bampair) {
            return "--bam";
          } else {
            return null;
          }
        }

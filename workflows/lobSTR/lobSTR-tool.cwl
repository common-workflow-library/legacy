#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: lobSTR
#"doap:homepage": http://melissagymrek.com/lobstr-code/index.html
#"doap:license": http://spdx.org/licenses/GPL-3.0
requirements:
- class: InlineJavascriptRequirement

inputs:

  maq:
    type: int
    default: 100
    inputBinding:
      prefix: --mapq
    doc: maximum allowed mapq score calculated as the sum of qualities at base mismatches
  minflank:
    type: int
    default: 8
    inputBinding:
      prefix: --minflank
    doc: Minimum length of flanking region to try to align
  bwaq:
    type: int?
    inputBinding:
      prefix: --bwaq
    doc: trim reads based on quality score. Same as the -q parameter in BWA
  reference:

    type: File
    #format: "edam:format_1929"
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
    - ${return self.location.replace(/(.*)ref\.fasta/, "$1chromsizes.tab");}
    - ${return self.location.replace(/(.*)ref\.fasta/, "$1mergedref.bed");}
    - ${return self.location.replace(/(.*)ref\.fasta/, "$1ref_map.tab");}

    doc: lobSTR's bwa reference files
  min-read-length:
    type: int
    default: 45
    inputBinding:
      prefix: --min-read-length
    doc: Don't process reads shorter than this
  bampair:
    type: boolean?
    inputBinding:
      prefix: --bampair
    doc: |
      Reads are in BAM format and are paired-end. NOTE: BAM file MUST be sorted or collated by read name.
  min-flank-allow-mismatch:
    type: int
    default: 30
    inputBinding:
      prefix: --min-flank-allow-mismatch
    doc: Don't allow mismatches if aligning flanking regions shorter than this
  rg-lib:
    type: string
    inputBinding:
      prefix: --rg-lib
    doc: Use this in the read group LB tag
  max-read-length:
    type: int
    default: 1024
    inputBinding:
      prefix: --min-read-length
    doc: Don't process reads longer than this
  r:
    type: float
    default: -1
    inputBinding:
      prefix: -e
    doc: fraction of missing alignments given 2% uniform base error rate. Ignored
      if --mismatch is set
  max-hits-quit-aln:
    type: int
    default: 1000
    inputBinding:
      prefix: --max-hits-quit-aln
    doc: Stop alignment search after this many hits found. Use -1 for no limit.
  files:
    type: File[]?
    #format: ["edam:format_1929", "edam:format_1930", "edam:format_2572"]
    inputBinding:
      prefix: --files
      position: 2
      itemSeparator: ','
    doc: List of files of single ends in fasta, fastq, or BAM files to align. Also
      use this parameter to specify paired-end BAM files to align.
  fft-window-step:
    type: int
    default: 4
    inputBinding:
      prefix: --fft-window-step
    doc: Step size of the sliding window
  extend:
    type: int
    default: 1000
    inputBinding:
      prefix: --extend
    doc: Number of bp the reference was extended when building the index. Must be
      same as --extend parameter used to run lobstr_index.py
  mismatch:
    type: int
    default: 1
    inputBinding:
      prefix: --mismatch
    doc: allowed edit distance in either flanking region
  max-diff-ref:
    type: int
    default: 50
    inputBinding:
      prefix: --max-diff-ref
    doc: Only report reads differing by at most this number of bp from the reference
      allele
  oldillumina:
    type: boolean?
    inputBinding:
      prefix: --oldillumina
    doc: |
      specifies that base pair quality scores are reported in the old Phred
      format (Illumina 1.3+, Illumina 1.5+) where quality scores are given as
      Phred + 64 rather than Phred + 33
  maxflank:
    type: int
    default: 100
    inputBinding:
      prefix: --maxflank
    doc: Length to trim flanking regions to before aligning
  e:
    type: int
    default: 1
    inputBinding:
      prefix: -e
    doc: Maximum number of gap extends allowed in each flanking region
  p2:

    type: File[]?
    #format: ["edam:format_1929", "edam:format_1930"]
    inputBinding:
      prefix: --p2
      position: 3
      itemSeparator: ','
    doc: list of files containing the second end of paired end reads in fasta or fastq
      format
  multi:
    type: boolean?
    inputBinding:
      prefix: --multi
    doc: Report reads mapping to multiple genomic locations
  p1:

    type: File[]?
    #format: ["edam:format_1929", "edam:format_1930"]
    inputBinding:
      prefix: --p1
      position: 2
      itemSeparator: ','
    doc: list of files containing the first end of paired end reads in fasta or fastq
      format
  entropy-threshold:
    type: float
    default: 0.45
    inputBinding:
      prefix: --entropy-threshold
    doc: Threshold score to call a window periodic
  g:
    type: int
    default: 1
    inputBinding:
      prefix: -g
    doc: Maximum number of gap opens allowed in each flanking region
  fft-window-size:
    type: int
    default: 16
    inputBinding:
      prefix: --fft-window-size
    doc: Size of sliding entropy window
  output_prefix:

    type: string
    inputBinding:
      prefix: --out
    doc: prefix for output files. will output prefix.aligned.bam and prefix.aligned.stats
  u:
    type: boolean?
    inputBinding:
      prefix: -u
    doc: only report reads different by an integer number of copy numbers from the
      reference allele
  rg-sample:
    type: string
    inputBinding:
      prefix: --rg-sample
    doc: Use this in the read group SM tag
outputs:
  bam:
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.aligned.bam')

  bam_stats:
    type: File
    outputBinding:
      glob: $(inputs['output_prefix'] + '.aligned.stats')

baseCommand: [lobSTR]
arguments:
- --verbose
- --noweb
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
doc: lobSTR is a tool for profiling Short Tandem Repeats (STRs) from high throughput
  sequencing data.


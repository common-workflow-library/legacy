#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      trimmomatic:
        specs: [ "https://identifiers.org/rrid/RRID:SCR_011848" ]
        version: [ "0.32", "0.35", "0.36" ]

requirements:
 ResourceRequirement:
   ramMin: 10240
   coresMin: 8
 SchemaDefRequirement:
   types:
    - $import: trimmomatic-end_mode.yaml
    - $import: trimmomatic-sliding_window.yaml
    - $import: trimmomatic-phred.yaml
    - $import: trimmomatic-illumina_clipping.yaml
    - $import: trimmomatic-max_info.yaml
 InlineJavascriptRequirement: {}
 ShellCommandRequirement: {}

# hints:
#  - $import: trimmomatic-docker.yml

inputs:
  phred:
    type: trimmomatic-phred.yaml#phred?
    inputBinding:
      prefix: -phred
      separate: false
      position: 4
    doc: |
      "33" or "64" specifies the base quality encoding. Default: 64

  tophred64:
    type: boolean?
    inputBinding:
      position: 12
      prefix: TOPHRED64
      separate: false
    doc: This (re)encodes the quality part of the FASTQ file to base 64.

  headcrop:
    type: int?
    inputBinding:
      position: 13
      prefix: 'HEADCROP:'
      separate: false
    doc: |
      Removes the specified number of bases, regardless of quality, from the
      beginning of the read.
      The number specified is the number of bases to keep, from the start of
      the read.

  tophred33:
    type: boolean?
    inputBinding:
      position: 12
      prefix: TOPHRED33
      separate: false
    doc: This (re)encodes the quality part of the FASTQ file to base 33.

  minlen:
    type: int?
    inputBinding:
      position: 100
      prefix: 'MINLEN:'
      separate: false
    doc: |
      This module removes reads that fall below the specified minimal length.
      If required, it should normally be after all other processing steps.
      Reads removed by this step will be counted and included in the "dropped
      reads" count presented in the trimmomatic summary.

  java_opts:
    type: string?
    inputBinding:
      position: 1
      shellQuote: false
    doc: |
      JVM arguments should be a quoted, space separated list
      (e.g. "-Xms128m -Xmx512m")

  leading:
    type: int?
    inputBinding:
      position: 14
      prefix: 'LEADING:'
      separate: false
    doc: |
      Remove low quality bases from the beginning. As long as a base has a
      value below this threshold the base is removed and the next base will be
      investigated.

  slidingwindow:
    type: trimmomatic-sliding_window.yaml#slidingWindow?
    inputBinding:
      position: 15
      valueFrom: |
        ${ if ( self ) {
             return "SLIDINGWINDOW:" + self.windowSize + ":"
               + self.requiredQuality;
           } else {
             return self;
           }
         }
    doc: |
      Perform a sliding window trimming, cutting once the average quality
      within the window falls below a threshold. By considering multiple
      bases, a single poor quality base will not cause the removal of high
      quality data later in the read.
      <windowSize> specifies the number of bases to average across
      <requiredQuality> specifies the average quality required

  illuminaClip:
    type: trimmomatic-illumina_clipping.yaml#illuminaClipping?
    inputBinding:
      valueFrom: |
        ${ if ( self ) { 
            var i="ILLUMINACLIP:" + inputs.illuminaClip.adapters.path + ":"
              + self.seedMismatches + ":" + self.palindromeClipThreshold + ":"
              + self.simpleClipThreshold;
            if ( self.minAdapterLength && !self.keepBothReads ) {
              return i + ":" + self.minAdapterLength;
            } else if ( self.minAdapterLength && self.keepBothReads ) {
              return i + ":" + self.minAdapterLength + ":" + self.keepBothReads;
            } else if ( self.keepBothReads && !self.minAdapterLength ) {
              return i + ":8:" + self.keepBothReads;
            } else {
              return i;
            }
          } else {
            return self;
          }
        }
      position: 11
    doc: Cut adapter and other illumina-specific sequences from the read.

  crop:
    type: int?
    inputBinding:
      position: 13
      prefix: 'CROP:'
      separate: false
    doc: |
      Removes bases regardless of quality from the end of the read, so that the
      read has maximally the specified length after this step has been
      performed. Steps performed after CROP might of course further shorten the
      read. The value is the number of bases to keep, from the start of the read.

  reads2:
    type: File?
    format: edam:format_1930  # fastq
    inputBinding:
      position: 6
    doc: FASTQ file of R2 reads in Paired End mode

  reads1:
    type: File
    format: edam:format_1930  # fastq
    inputBinding:
      position: 5
    doc: FASTQ file of reads (R1 reads in Paired End mode)

  avgqual:
    type: int?
    inputBinding:
      position: 101
      prefix: 'AVGQUAL:'
      separate: false
    doc: |
      Drop the read if the average quality is below the specified level

  trailing:
    type: int?
    inputBinding:
      position: 14
      prefix: 'TRAILING:'
      separate: false
    doc: |
      Remove low quality bases from the end. As long as a base has a value
      below this threshold the base is removed and the next base (which as
      trimmomatic is starting from the 3' prime end would be base preceding
      the just removed base) will be investigated. This approach can be used
      removing the special Illumina "low quality segment" regions (which are
      marked with quality score of 2), but we recommend Sliding Window or
      MaxInfo instead

  maxinfo:
    type: trimmomatic-max_info.yaml#maxinfo?
    inputBinding:
      position: 15
      valueFrom: |
        ${ if ( self ) {
             return "MAXINFO:" + self.targetLength + ":" + self.strictness;
           } else {
             return self;
           }
         }
    doc: |
      Performs an adaptive quality trim, balancing the benefits of retaining
      longer reads against the costs of retaining bases with errors.
      <targetLength>: This specifies the read length which is likely to allow
      the location of the read within the target sequence to be determined.
      <strictness>: This value, which should be set between 0 and 1, specifies
      the balance between preserving as much read length as possible vs.
      removal of incorrect bases. A low value of this parameter (<0.2) favours
      longer reads, while a high value (>0.8) favours read correctness.

  end_mode:
    type: trimmomatic-end_mode.yaml#end_mode
    inputBinding:
      position: 3
    doc: |
      Single End (SE) or Paired End (PE) mode

outputs:
  reads1_trimmed:
    type: File
    format: edam:format_1930  # fastq
    outputBinding:
      glob: $(inputs.reads1.nameroot).trimmed.fastq

  output_log:
    type: File
    outputBinding:
      glob: trim.log
    label: Trimmomatic log
    doc: |
      log of all read trimmings, indicating the following details:
        the read name
        the surviving sequence length
        the location of the first surviving base, aka. the amount trimmed from the start
        the location of the last surviving base in the original read
        the amount trimmed from the end

  reads1_trimmed_unpaired:
    type: File?
    format: edam:format_1930  # fastq
    outputBinding:
      glob: $(inputs.reads1.nameroot).trimmed.unpaired.fastq

  reads2_trimmed_paired:
    type: File?
    format: edam:format_1930  # fastq
    outputBinding:
      glob: |
        ${ if (inputs.reads2 ) {
             return inputs.reads2.nameroot + '.trimmed.fastq';
           } else {
             return null;
           }
         }

  reads2_trimmed_unpaired:
    type: File?
    format: edam:format_1930  # fastq
    outputBinding:
      glob: |
        ${ if (inputs.reads2 ) {
             return inputs.reads2.nameroot + '.trimmed.unpaired.fastq';
           } else {
             return null;
           }
         }

baseCommand: [ java, org.usadellab.trimmomatic.Trimmomatic ]

arguments:
- valueFrom: trim.log
  prefix: -trimlog 
  position: 4
- valueFrom: $(runtime.cores)
  position: 4
  prefix: -threads
- valueFrom: $(inputs.reads1.nameroot).trimmed.fastq
  position: 7
- valueFrom: |
    ${
      if (inputs.end_mode == "PE" && inputs.reads2) {
        return inputs.reads1.nameroot + '.trimmed.unpaired.fastq';
      } else {
        return null;
      }
    }
  position: 8
- valueFrom: |
    ${
      if (inputs.end_mode == "PE" && inputs.reads2) {
        return inputs.reads2.nameroot + '.trimmed.fastq';
      } else {
        return null;
      }
    }
  position: 9
- valueFrom: |
    ${
      if (inputs.end_mode == "PE" && inputs.reads2) {
        return inputs.reads2.nameroot + '.trimmed.unpaired.fastq';
      } else {
        return null;
      }
    }
  position: 10

doc: |
  Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop
  Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
  depending on the library preparation and downstream application.
  There are two major modes of the program: Paired end mode and Single end mode. The
  paired end mode will maintain correspondence of read pairs and also use the additional
  information contained in paired reads to better find adapter or PCR primer fragments
  introduced by the library preparation process.
  Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores,
  depending on the Illumina pipeline used).

$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"

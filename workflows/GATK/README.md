# GATK h3abionet pipeline docs
# Overview

A [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS) germline workflow designed to work with GATK 3.5  (Van der Auwera et al., 2013).

# Dependencies
The pipeline rely on the following requirements.

## System

## Tools
### FastQC
FastQC is used as an initial QC step where the input files are checked for usual metrics such as:
  - Read length
  - Reads distribution
  - GC content
  - ...

### Trimmomatic
Trimmomatic is the entry point of the pipeline, it is used to cleanup the reads in the input fastq files from any sequencing adaptors.

### BWA
BWA is used to align the reads from the the input fastq files (paired-ends). We use specifically `bwa mem` as recommended by the [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS). BWA produces a SAM file containing the aligned reads against the human reference genome (hg19, GATK bundle build 2.8).

As GATK tools downstream requires properly formatted Read Group information. We add by default 'toy' Read Group information while processing the alignment to the output SAM file. we specifically use the flag `-R '@RG\tID:foo\tSM:bar\tLB:library1'`.

### SAMtools
[SAMtools](http://www.htslib.org) (Li et al., 2009) are used few times in the pipeline:
  1. Convert BWA's output from a SAM format to a BAM format
  2. Sort the reads in the generated BAM file in step 1 (above)
  3. Indexing the BAM file for the following tools to use

### Picard
[Picard tools](https://broadinstitute.github.io/picard/) are used to mark duplicated reads in the aligned and sorted BAM file, making thus the files lighter and less prone to errors in the downstream steps of the pipeline.

### GATK
(DePristo et al., 2011)
### SNPEff
### BAMStat


## Reference Files

# Running

# TODO
- [ ] streaming between steps to improve throughput
   - [ ] need to stream between samtools/bwa steps for performance
- [ ] properly configuring/exposing threads per step to maximize speed of the steps
- [ ] remove the dictionary creation steps
- [ ] clearly document the reference files needed
- [ ] dealing with the CWL input file directory being not writable (this prevents GATK from making VCF index files on the fly). The error looks like this: `WARN  09:30:55,672 RMDTrackBuilder - Unable to write to /var/lib/cwl/stg24cab595-694a-4cac-b90a-abe40abdeaa2/Mills_and_1000G_gold_st.vcf.idx for the index file, creating index in memory only`
- [ ] adding SV and CNV calling tools to the workflow
 	 - [ ] CNV: [CNVnator](http://sv.gersteinlab.org/)
 	 - [ ] SV: [Delly2](https://github.com/tobiasrausch/delly)
- [ ] adding Oncotator and Annovar for annotation (have SNPEff)
- [ ] report error back to CWLTools where schemas are loaded on each run (should be cached) and if there's a failure to retrieve an 'utf8' decoding error is reported.
- [ ] output cleanup, we don't need to save all the files the workflow currently does
- [ ] zip file for DepthOfCoverage

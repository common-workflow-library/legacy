# Overview

A GATK best-practices germline workflow designed to work with GATK 3.5.

# Dependencies

## System

## Reference Files

# Running

We use the CWL reference implementation.

    cwltool --cachedir <path_to_tmp_dir> --debug --non-strict --tmpdir-prefix  GATK-complete-Workflow-h3abionet.cwl GATK-complete-Workflow-copy-h3abionet.json

# TODO

- [ ] streaming between steps to improve throughput, currently not supported via the reference implementation
   - [ ] need to stream between samtools/bwa steps for performance
- [X] properly configuring/exposing threads per step to maximize speed of the steps
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
- [ ] need to stream between samtools/bwa steps for performance
- [ ] add scatter/gather based on chr to the workflow
- [ ] DepthOfCoverage need to have the input from view and sort steps in the same directory (e.g. the bam and bai need to be mounted in the same directory). CWL puts each input on it's own path and these are read only.  Michael is going to help us work around this.  In the mean time, DepthOfCoverage is commented out

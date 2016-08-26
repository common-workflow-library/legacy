# Overview

A GATK best-practices germline workflow designed to work with GATK 3.5.

# Dependencies

## System

## Reference Files

# Running

# TODO

* streaming between steps to improve throughput
* properly configuring/exposing threads per step to maximize speed of the steps
* remove the dictionary creation steps
* clearly document the reference files needed
* dealing with the CWL input file directory being not writable (this prevents GATK from making VCF index files on the fly)
* adding SV and CNV calling tools to the workflow
  * CNV: [CNVnator](http://sv.gersteinlab.org/)
  * SV: [Delly2](https://github.com/tobiasrausch/delly)
* adding Oncotator and Annovar for annotation (have SNPEff)

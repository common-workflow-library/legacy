Genome Analysis ToolKit (GATK) best Practices workflows provide step-by-step recommendations for performing variant discovery analysis in high-throughput sequencing (HTS) data. This directory is example implementation of GATK varaint calling workflow in CWL. It has been designed to deal with Forward and reverse read files (provided in test-data). 
Following pre-processinf steps are required: 

1. To create index of hg19.fa run : bwa index -a bwtsw hg19.fa. hg19.fa is provided in test-data. 
2. To create hg19.fa.fai, run: samtools faidx hg19.fa

Please note that wrappers for both tools (bwa index and samtools faidx) are present in CWL/workflows/tools repo and can be used for the above purpose if bwa and samtools are not installed on the local system. 
The purpose of keeping these two steps for preprocessing instead of including in workflow is that these files can be reused again and again as hg19.fa does not change if input files change. 

To run locally, you may install the reference implementation of CWL:

1. Install the reference implementation of cwl-runner with  "pip install cwl-runner".
2. Download "test-data" from "https://www.dropbox.com/sh/0bpamtulxgdcf11/AACM8bpbT87BKudkEv8VGOXwa?dl=0" and save the directory in working directory. This test data contains 3 known variant files, sample forward read and reverse read file for NA12878 and compressed hg19.fa.zip and enough to run the given workflow. 
3. Create a temporary directory named "tmpdir" in current working directory for further processing.
4. Run the demo available with the following command: 
   cwltool --tmpdir-prefix=$(pwd)/tmpdir --tmp-outdir-prefix=$(pwd)/tmpdir ./GATK-workflow.cwl ./GATK-workflow.json


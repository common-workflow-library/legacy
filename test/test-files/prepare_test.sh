#!/usr/bin/env bash

# Script has to be run in current directory ./test/test-files

#download test files
wget -q -O - ftp://hgdownload.cse.ucsc.edu/goldenPath/dm3/chromosomes/chr4.fa.gz|gunzip >dm3_chr4.fa
wget -q -O - https://scidap.com/tests/dm3.gtf.bz2 | bunzip2 > dm3_chr4.gtf
wget -q -O - https://scidap.com/tests/SRR1031972.fastq.bz2 |bunzip2 > SRR1031972.fastq
mkdir -p dm3

#generate required files. Reference genome. Align reads. Create index.
cd ..
cwltool ../tools/STAR.cwl ./STAR-genomeGenerate-job.json 
cwltool ../tools/STAR.cwl ./STAR-alignReads-job.json
cwltool ../tools/samtools-index.cwl ./samtools-index-job.json 
mv SRR1031972.Aligned.sortedByCoord.out.bam.bai ./test-files/
cd test-files

#link the job file into current directory. and run the pipeline
ln -s ../../workflows/scidap/bam-genomecov-bigwig-job.json  ./bam-genomecov-bigwig-job.json 
cwltool --tmpdir-prefix $(pwd) --tmp-outdir-prefix $(pwd) --debug ../../workflows/scidap/bam-genomecov-bigwig.cwl ./bam-genomecov-bigwig-job.json 
rm ./bam-genomecov-bigwig-job.json

Sample CWL Command Line Tools.

# Testing CWLs

Test directory includes:
* dm3_chr4.fa - Chromosome 4 of Drosophila genome
* dm3_chr4.gtf - Chromosome 4 RefSeq annotation file
* SRR1031972.fastq - The reduced raw reads file ( reads from Chromosome 4 only, RNA-Seq data)

To test tools they have to be executed in particular order to produce input for the next one. Make ```./workflow/tools``` your current working directory 
and you will run command like this ```cwtool --basedir ./ ./TOOL.cwl ./jobs/TOOL-job.json```.

## Tools

Indexing genome
---------------

The first step is to indexing our genome. I'm going to use *STAR genereteGenome*. To do that 
I run ```cwltool --basedir ./ ./STAR.cwl ./jobs/STAR-job-index.json```. dm3 genome will be place into ./test-files/dm3 directory.

Output from cwltool:

```
/usr/local/bin/cwltool 1.0.20151125221324
[job 4532697104] ./workflows/tools$ docker run -i --volume=./workflows/tools/test-files/dm3_chr4.gtf:/tmp/job257401672_test-files/dm3_chr4.gtf:ro --volume=./workflows/tools/test-files/dm3_chr4.fa:/tmp/job257401672_test-files/dm3_chr4.fa:ro --volume=./workflows/tools:/tmp/job_output:rw --volume=/var/folders/hx/3qsmpl9s50zdmn49jb03l_tw0000gn/T/tmpWo01Wd:/tmp/job_tmp:rw --workdir=/tmp/job_output --read-only=true --user=1000 --rm --env=TMPDIR=/tmp/job_tmp --env=PATH=/usr/local/bin/:/usr/bin:/bin scidap/star:v2.5.0a STAR --genomeDir ./test-files/dm3/ --genomeFastaFiles /tmp/job257401672_test-files/dm3_chr4.fa --outBAMcompression 10 --outSAMmode Full --outSAMtype BAM SortedByCoordinate --outStd Log --runMode genomeGenerate --runThreadN 4 --sjdbGTFfile /tmp/job257401672_test-files/dm3_chr4.gtf --sjdbOverhang 100
Nov 26 17:57:46 ..... Started STAR run
Nov 26 17:57:46 ... Starting to generate Genome files
Nov 26 17:57:46 ... starting to sort  Suffix Array. This may take a long time...
Nov 26 17:57:46 ... sorting Suffix Array chunks and saving them to disk...
Nov 26 17:57:47 ... loading chunks from disk, packing SA...
Nov 26 17:57:47 ... Finished generating suffix array
Nov 26 17:57:47 ... Generating Suffix Array index
Nov 26 17:57:50 ... Completed Suffix Array index
Nov 26 17:57:50 ..... Processing annotations GTF
Nov 26 17:57:50 ..... Inserting junctions into the genome indices
Nov 26 17:57:58 ... writing Genome to disk ...
Nov 26 17:57:58 ... writing Suffix Array to disk ...
Nov 26 17:57:58 ... writing SAindex to disk
Nov 26 17:58:35 ..... Finished successfully
Final process status is success
{
    "indices": {
        "path": "./test-files/dm3//Genome",
        "size": 1753563,
        "secondaryFiles": [
            {
                "path": "./test-files/dm3//SA",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//SAindex",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//chrNameLength.txt",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//chrLength.txt",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//chrStart.txt",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//geneInfo.tab",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//sjdbList.fromGTF.out.tab",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//chrName.txt",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//exonGeTrInfo.tab",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//genomeParameters.txt",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//sjdbList.out.tab",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//exonInfo.tab",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//sjdbInfo.txt",
                "class": "File"
            },
            {
                "path": "./test-files/dm3//transcriptInfo.tab",
                "class": "File"
            }
        ],
        "class": "File",
        "checksum": "sha1$761906d19ceb10a0e2677afdfb756c4f1ca925a1"
    },
    "aligned": null,
    "mappingstats": null
}%
```

Reads alignment
---------------

To align reads run ```cwltool --basedir ./ ./STAR.cwl ./jobs/STAR-job-rna.json```

```
/usr/local/bin/cwltool 1.0.20151126171959
[job 4518994960] ./workflows/tools$ docker run -i --volume=./workflows/tools/test-files/dm3/:/tmp/job958208261_test-files/dm3/:ro --volume=./workflows/tools/test-files/SRR1031972.fastq:/tmp/job958208261_test-files/SRR1031972.fastq:ro --volume=./workflows/tools:/tmp/job_output:rw --volume=/var/folders/hx/3qsmpl9s50zdmn49jb03l_tw0000gn/T/tmptuCn_d:/tmp/job_tmp:rw --workdir=/tmp/job_output --read-only=true --user=1000 --rm --env=TMPDIR=/tmp/job_tmp --env=PATH=/usr/local/bin/:/usr/bin:/bin scidap/star:v2.5.0a STAR --genomeDir /tmp/job958208261_test-files/dm3/ --outBAMcompression 10 --outFileNamePrefix ./test-files/SRR1031972. --outSAMmode Full --outSAMtype BAM SortedByCoordinate --outStd Log --readFilesIn /tmp/job958208261_test-files/SRR1031972.fastq --runMode alignReads --runThreadN 4
Nov 26 19:27:58 ..... Started STAR run
Nov 26 19:27:58 ..... Loading genome
Nov 26 19:28:08 ..... Started mapping
Nov 26 19:29:34 ..... Started sorting BAM
Nov 26 19:29:36 ..... Finished successfully
Final process status is success
{
    "indices": null, 
    "aligned": {
        "path": "./workflows/tools/./test-files/SRR1031972.Aligned.sortedByCoord.out.bam", 
        "size": 22139153, 
        "secondaryFiles": [
            {
                "path": "./test-files/SRR1031972.Log.final.out", 
                "class": "File"
            }, 
            {
                "path": "./test-files/SRR1031972.SJ.out.tab", 
                "class": "File"
            }, 
            {
                "path": "./test-files/SRR1031972.Log.out", 
                "class": "File"
            }
        ], 
        "class": "File", 
        "checksum": "sha1$ba38fcd1f238553d244f339d1147cd591324e207"
    }, 
    "mappingstats": "[{\"Started job on \":\"Nov 26 19:27:58\"},{\"Started mapping on \":\"Nov 26 19:28:08\"},{\"Finished on \":\"Nov 26 19:29:36\"},{\"Mapping speed, Million of reads per hour \":\"8.55\"},{\"Number of input reads \":\"209081\"},{\"Average input read length \":\"40\"},{\"Uniquely mapped reads number \":\"64313\"},{\"Uniquely mapped reads % \":\"30.76%\"},{\"Average mapped length \":\"38.19\"},{\"Number of splices: Total \":\"14213\"},{\"Number of splices: Annotated (sjdb) \":\"1640\"},{\"Number of splices: GT/AG \":\"12072\"},{\"Number of splices: GC/AG \":\"299\"},{\"Number of splices: AT/AC \":\"3\"},{\"Number of splices: Non-canonical \":\"1839\"},{\"Mismatch rate per base, % \":\"1.86%\"},{\"Deletion rate per base \":\"0.00%\"},{\"Deletion average length \":\"1.25\"},{\"Insertion rate per base \":\"0.00%\"},{\"Insertion average length \":\"1.03\"},{\"Number of reads mapped to multiple loci \":\"144768\"},{\"% of reads mapped to multiple loci \":\"69.24%\"},{\"Number of reads mapped to too many loci \":\"0\"},{\"% of reads mapped to too many loci \":\"0.00%\"},{\"% of reads unmapped: too many mismatches \":\"0.00%\"},{\"% of reads unmapped: too short \":\"0.00%\"},{\"% of reads unmapped: other \":\"0.00%\"}]"
}%                               
```

Indexing .bam file
------------------

To index the .bam file ```cwltool --basedir ./ --outdir ./test-files ./samtools-index.cwl ./jobs/samtools-index-job.json```

Result:
```json
{
    "sorted": {
        "path": "././test-files/SRR1031972.Aligned.sortedByCoord.out.bam.bai", 
        "size": 40528, 
        "class": "File", 
        "checksum": "sha1$83738ffada23f654ba1f471973c7dccceb14cffc"
    }
}
```

Genome coverage
---------------

To create a genome coverage file .bedGraph ```cwltool --basedir ./ ./bedtools-genomecov.cwl ./jobs/bedtools-genomecov-job.json```

Result:
```json
{
    "genomecoverage": {
        "path": "./workflows/tools/./test-files/SRR1031972.bedGraph", 
        "size": 1423143, 
        "class": "File", 
        "checksum": "sha1$dd87be96fc201734c2e5017f86e056b4bb0b2b3f"
    }
}      
```

# CWL Tools & Workflows

![CWL Logo](https://cloud.githubusercontent.com/assets/1330696/11716362/045022ce-9f52-11e5-827d-504773212449.jpg =100x100)


[![Build Status](https://travis-ci.org/common-workflow-language/workflows.svg?branch=master)](https://travis-ci.org/common-workflow-language/workflows)


"CWL Tools & Workflows" is the community's best practices for CWL tools and workflows description. 

To submit a new workflow or tool description please make a pull request against this repository.

Any contribution submitted for inclusion in this repository shall be under the
terms and conditions of the Apache License version 2.0 (see LICENSE.txt).

## Pull request

Please consider these recommendations before pull request:

* One file per tool or per logical part of the tool like with subcommands
* CWL’s filename is the same as tool’s name (```STAR.cwl```)  if possible and no conflicts
* If a tool has subcommands, like samtools or bwa use tool name dash subcommand name. For example, ```bwa-mem.cwl``` or ```samtools-index.cwl```
* Place CWL files into tools directory
* For each tool create job file and test file and place them into test directory
 * Filename for the job file has to have the same basename plus *-job.json* (bwa-mem-job.json, samtools-index-job.json)
 * Filename for the test file has to have the same basename plus *-test.yaml* (bwa-mem-test.yaml, samtools-index-test.yaml)
* Please, use docker for the tool your are describing
  * If the tool has subcommands and you will have the same code in all files. Please, move it into separate file (```samtools-docker.cwl```) and include it
  * If you are a maintainer of the image, it is good to provide content of Dockerfile in class DockerRequirement (dockerFile)

Incomplete descriptions are welcome as long as they are usable.
General encouragement to share early & often

## Testing CWLs

Test directory includes:
* dm3_chr4.fa - Chromosome 4 of Drosophila genome
* dm3_chr4.gtf - Chromosome 4 RefSeq annotation file
* SRR1031972.fastq - The reduced raw reads file ( reads from Chromosome 4 only, RNA-Seq data)

To test tools run ```test/test.sh``` from the repository root.


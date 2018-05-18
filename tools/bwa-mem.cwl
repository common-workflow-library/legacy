#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: bwa-docker.yml
- class: InlineJavascriptRequirement

inputs:
  minimum_seed_length:
    type: int?
    inputBinding:
      position: 1
      prefix: -k
    doc: -k INT        minimum seed length [19]
  reference:
    type: File
    secondaryFiles:
        - .amb
        - .ann
        - .bwt
        - .pac
        - .sa
    inputBinding:
      position: 2

  output_filename: string
  reads:
    type:
      type: array
      items: File
    inputBinding:
      position: 3

  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: -t
    doc: -t INT        number of threads [1]
  min_std_max_min:
    type: int[]?
    inputBinding:
      position: 1
      prefix: -I
      itemSeparator: ','
stdout: $(inputs.output_filename)

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

baseCommand:
- bwa
- mem

doc: |
  Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]

  Algorithm options:
         -w INT        band width for banded alignment [100]
         -d INT        off-diagonal X-dropoff [100]
         -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
         -y INT        seed occurrence for the 3rd round seeding [20]
         -c INT        skip seeds with more than INT occurrences [500]
         -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
         -W INT        discard a chain if seeded bases shorter than INT [0]
         -m INT        perform at most INT rounds of mate rescues for each read [50]
         -S            skip mate rescue
         -P            skip pairing; mate rescue performed unless -S also in use
         -e            discard full-length exact matches

  Scoring options:

         -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
         -B INT        penalty for a mismatch [4]
         -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
         -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
         -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
         -U INT        penalty for an unpaired read pair [17]

         -x STR        read type. Setting -x changes multiple parameters unless overriden [null]
                       pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
                       ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
                       intractg: -B9 -O16 -L5  (intra-species contigs to ref)

  Input/output options:

         -p            smart pairing (ignoring in2.fq)
         -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
         -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
         -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)

         -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
         -T INT        minimum score to output [30]
         -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
         -a            output all alignments for SE or unpaired PE
         -C            append FASTA/FASTQ comment to SAM output
         -V            output the reference FASTA header in the XR tag
         -Y            use soft clipping for supplementary alignments
         -M            mark shorter split hits as secondary

         -I FLOAT[,FLOAT[,INT[,INT]]]
                       specify the mean, standard deviation (10% of the mean if absent), max
                       (4 sigma from the mean if absent) and min of the insert size distribution.
                       FR orientation only. [inferred]

  Note: Please read the man page for detailed description of the command line and options.


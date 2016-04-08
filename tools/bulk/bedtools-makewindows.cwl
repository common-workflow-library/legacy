#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: bedtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#g'
    type:
      - 'null'
      - boolean
    description: |
      <genome>
      Genome file size (see notes below).
      Windows will be created for each chromosome in the file.
    inputBinding:
      position: 1
      prefix: '-g'
  - id: '#b'
    type:
      - 'null'
      - boolean
    description: |
      <bed>
      BED file (with chrom,start,end fields).
      Windows will be created for each interval in the file.
      Windows Output Options:
    inputBinding:
      position: 1
      prefix: '-b'
  - id: '#w'
    type:
      - 'null'
      - boolean
    description: |
      <window_size>
      Divide each input interval (either a chromosome or a BED interval)
      to fixed-sized windows (i.e. same number of nucleotide in each window).
      Can be combined with -s <step_size>
    inputBinding:
      position: 1
      prefix: '-w'
  - id: '#s'
    type:
      - 'null'
      - boolean
    description: |
      <step_size>
      Step size: i.e., how many base pairs to step before
      creating a new window. Used to create "sliding" windows.
      - Defaults to window size (non-sliding windows).
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#n'
    type:
      - 'null'
      - boolean
    description: |
      <number_of_windows>
      Divide each input interval (either a chromosome or a BED interval)
      to fixed number of windows (i.e. same number of windows, with
      varying window sizes).
      -reverse
      Reverse numbering of windows in the output, i.e. report
      windows in decreasing order
      ID Naming Options:
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#i'
    type:
      - 'null'
      - boolean
    description: |
      src|winnum|srcwinnum
      The default output is 3 columns: chrom, start, end .
      With this option, a name column will be added.
      "-i src" - use the source interval's name.
      "-i winnum" - use the window number as the ID (e.g. 1,2,3,4...).
      "-i srcwinnum" - use the source interval's name with the window number.
      See below for usage examples.
    inputBinding:
      position: 1
      prefix: '-i'
outputs:
  - id: '#stdoutfile'
    type: File
    outputBinding:
      glob:
        engine: 'cwl:JsonPointer'
        script: /job/stdoutfile
stdout:
  engine: 'cwl:JsonPointer'
  script: /job/stdoutfile
baseCommand:
  - bedtools
  - makewindows
description: |
  Tool: bedtools makewindows
  Version: v2.25.0
  Summary: Makes adjacent or sliding windows across a genome or BED file.

  Usage: bedtools makewindows [OPTIONS] [-g <genome> OR -b <bed>]
   [ -w <window_size> OR -n <number of windows> ]

  Input Options: 
  	-g <genome>
  		Genome file size (see notes below).
  		Windows will be created for each chromosome in the file.

  	-b <bed>
  		BED file (with chrom,start,end fields).
  		Windows will be created for each interval in the file.

  Windows Output Options: 
  	-w <window_size>
  		Divide each input interval (either a chromosome or a BED interval)
  		to fixed-sized windows (i.e. same number of nucleotide in each window).
  		Can be combined with -s <step_size>

  	-s <step_size>
  		Step size: i.e., how many base pairs to step before
  		creating a new window. Used to create "sliding" windows.
  		- Defaults to window size (non-sliding windows).

  	-n <number_of_windows>
  		Divide each input interval (either a chromosome or a BED interval)
  		to fixed number of windows (i.e. same number of windows, with
  		varying window sizes).

  	-reverse
  		 Reverse numbering of windows in the output, i.e. report 
  		 windows in decreasing order

  ID Naming Options: 
  	-i src|winnum|srcwinnum
  		The default output is 3 columns: chrom, start, end .
  		With this option, a name column will be added.
  		 "-i src" - use the source interval's name.
  		 "-i winnum" - use the window number as the ID (e.g. 1,2,3,4...).
  		 "-i srcwinnum" - use the source interval's name with the window number.
  		See below for usage examples.

  Notes: 
  	(1) The genome file should tab delimited and structured as follows:
  	 <chromName><TAB><chromSize>

  	For example, Human (hg19):
  	chr1	249250621
  	chr2	243199373
  	...
  	chr18_gl000207_random	4262

  Tips: 
  	One can use the UCSC Genome Browser's MySQL database to extract
  	chromosome sizes. For example, H. sapiens:

  	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
  	"select chrom, size from hg19.chromInfo" > hg19.genome

  Examples: 
   # Divide the human genome into windows of 1MB:
   $ bedtools makewindows -g hg19.txt -w 1000000
   chr1 0 1000000
   chr1 1000000 2000000
   chr1 2000000 3000000
   chr1 3000000 4000000
   chr1 4000000 5000000
   ...

   # Divide the human genome into sliding (=overlapping) windows of 1MB, with 500KB overlap:
   $ bedtools makewindows -g hg19.txt -w 1000000 -s 500000
   chr1 0 1000000
   chr1 500000 1500000
   chr1 1000000 2000000
   chr1 1500000 2500000
   chr1 2000000 3000000
   ...

   # Divide each chromosome in human genome to 1000 windows of equal size:
   $ bedtools makewindows -g hg19.txt -n 1000
   chr1 0 249251
   chr1 249251 498502
   chr1 498502 747753
   chr1 747753 997004
   chr1 997004 1246255
   ...

   # Divide each interval in the given BED file into 10 equal-sized windows:
   $ cat input.bed
   chr5 60000 70000
   chr5 73000 90000
   chr5 100000 101000
   $ bedtools makewindows -b input.bed -n 10
   chr5 60000 61000
   chr5 61000 62000
   chr5 62000 63000
   chr5 63000 64000
   chr5 64000 65000
   ...

   # Add a name column, based on the window number: 
   $ cat input.bed
   chr5  60000  70000 AAA
   chr5  73000  90000 BBB
   chr5 100000 101000 CCC
   $ bedtools makewindows -b input.bed -n 3 -i winnum
   chr5        60000   63334   1
   chr5        63334   66668   2
   chr5        66668   70000   3
   chr5        73000   78667   1
   chr5        78667   84334   2
   chr5        84334   90000   3
   chr5        100000  100334  1
   chr5        100334  100668  2
   chr5        100668  101000  3
   ...

   # Reverse window numbers: 
   $ cat input.bed
   chr5  60000  70000 AAA
   chr5  73000  90000 BBB
   chr5 100000 101000 CCC
   $ bedtools makewindows -b input.bed -n 3 -i winnum -reverse
   chr5        60000   63334   3
   chr5        63334   66668   2
   chr5        66668   70000   1
   chr5        73000   78667   3
   chr5        78667   84334   2
   chr5        84334   90000   1
   chr5        100000  100334  3
   chr5        100334  100668  2
   chr5        100668  101000  1
   ...

   # Add a name column, based on the source ID + window number: 
   $ cat input.bed
   chr5  60000  70000 AAA
   chr5  73000  90000 BBB
   chr5 100000 101000 CCC
   $ bedtools makewindows -b input.bed -n 3 -i srcwinnum
   chr5        60000   63334   AAA_1
   chr5        63334   66668   AAA_2
   chr5        66668   70000   AAA_3
   chr5        73000   78667   BBB_1
   chr5        78667   84334   BBB_2
   chr5        84334   90000   BBB_3
   chr5        100000  100334  CCC_1
   chr5        100334  100668  CCC_2
   chr5        100668  101000  CCC_3
   ...


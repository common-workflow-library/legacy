#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  #dockerImageId: scidap/bowtie:v1.1.2 #not yet ready
  dockerPull: scidap/bowtie:v1.1.2
  dockerFile: >
    $import: bowtie-Dockerfile
inputs:
  Q1:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --Q1
    doc: |
      --Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
  fullref:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --fullref
    doc: "write entire ref name (default: only up to 1st space)"
  verbose:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --verbose
    doc: verbose output (for debugging)
  sam:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -S
    doc: |
      --sam           write hits in SAM format
  col-cseq:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --col-cseq
    doc: print aligned colorspace seqs as colors, not decoded bases
  fr:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --fr
    doc: |
      --rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
  quiet:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --quiet
    doc: print nothing but the alignments
  pairtries:
    type: int?
    inputBinding:
      position: 1
      prefix: --pairtries
    doc: |
      <int>  max # attempts to find mate for anchor hit (default: 100)
  seed:
    type: int?
    inputBinding:
      position: 1
      prefix: --seed
    doc: |
      <int>       seed for random number generator
  o:
    type: int?
    inputBinding:
      position: 1
      prefix: -o
    doc: |
      --offrate <int> override offrate of index; must be >= index's offrate
  snpphred:
    type: int?
    inputBinding:
      position: 1
      prefix: --snpphred
    doc: |
      <int>   Phred penalty for SNP when decoding colorspace (def: 30)
      or
  sam-nohead:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --sam-nohead
    doc: supppress header lines (starting with @) for SAM output
  col-cqual:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --col-cqual
    doc: print original colorspace quals, not decoded quals
  best:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --best
    doc: hits guaranteed best stratum; ties broken by quality
  sam-RG:
    type: string?
    inputBinding:
      position: 1
      prefix: --sam-RG
    doc: |
      <text>    add <text> (usually "lab=value") to @RG line of SAM header
      Performance:
  col-keepends:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --col-keepends
    doc: keep nucleotides at extreme ends of decoded alignment
  filename:
    type: string
    inputBinding:
      position: 11

  '3':
    type: int?
    inputBinding:
      position: 1
      prefix: '-3'
    doc: |
      --trim3 <int>   trim <int> bases from 3' (right) end of reads
  un:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --un
    doc: |
      <fname>       write unaligned reads/pairs to file(s) <fname>
  '5':
    type: int?
    inputBinding:
      position: 1
      prefix: '-5'
    doc: |
      --trim5 <int>   trim <int> bases from 5' (left) end of reads
  maxbts:
    type: int?
    inputBinding:
      position: 1
      prefix: --maxbts
    doc: |
      <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
  suppress:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --suppress
    doc: |
      <cols>  suppresses given columns (comma-delim'ed) in default output
  sam-nosq:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --sam-nosq
    doc: supppress @SQ header lines for SAM output
  shmem:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --shmem
    doc: use shared mem for index; many bowties can share
  ebwt:
    type: string
    inputBinding:
      position: 8

    doc: The basename of the index to be searched. The basename is the name of any
      of the index files up to but not including the final .1.ebwt / .rev.1.ebwt /
      etc. bowtie looks for the specified index first in the current directory, then
      in the indexes subdirectory under the directory where the bowtie executable
      is located, then looks in the directory specified in the BOWTIE_INDEXES environment
      variable.
  mm:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --mm
    doc: use memory-mapped I/O for index; many bowties can share
  C:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -C
    doc: reads and index are in colorspace
  B:
    type: int?
    inputBinding:
      position: 1
      prefix: -B
    doc: |
      --offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
  y:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -y
    doc: |
      --tryhard       try hard to find valid alignments, at the expense of speed'
  I:
    type: int?
    inputBinding:
      position: 1
      prefix: -I
    doc: |
      --minins <int>  minimum insert size for paired-end alignment (default: 0)
  max:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --max
    doc: |
      <fname>      write reads/pairs over -m limit to file(s) <fname>
  M:
    type: int?
    inputBinding:
      position: 1
      prefix: -M
    doc: |
      <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
  strata:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --strata
    doc: "hits in sub-optimal strata aren't reported (requires --best)"
  al:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --al
    doc: |
      <fname>       write aligned reads/pairs to file(s) <fname>
  Q:
    type: File?
    inputBinding:
      position: 1
      prefix: -Q
    doc: |
      --quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f
  phred64-quals:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --phred64-quals
    doc: input quals are Phred+64 (same as --solexa1.3-quals)
  mapq:
    type: int?
    inputBinding:
      position: 1
      prefix: --mapq
    doc: |
      <int>       default mapping quality (MAPQ) to print for SAM alignments
  X:
    type: int?
    inputBinding:
      position: 1
      prefix: -X
    doc: |
      --maxins <int>  maximum insert size for paired-end alignment (default: 250)
  solexa1.3-quals:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --solexa1.3-quals
    doc: |
      input quals are from GA Pipeline ver. >= 1.3
  chunkmbs:
    type: int?
    inputBinding:
      position: 1
      prefix: --chunkmbs
    doc: |
      <int>   max megabytes of RAM for best-first search frames (def: 64)
  a:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -a
    doc: |
      --all           report all alignments per read (much slower than low -k)
  c:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -c
    doc: |
      query sequences given on cmd line (as <mates>, <singles>)
  e:
    type: int?
    inputBinding:
      position: 1
      prefix: -e
    doc: |
      --maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
  large-index:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --large-index
    doc: |
      force usage of a ''large'' index, even if a small one is present
  f:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -f
    doc: "query input files are (multi-)FASTA .fa/.mfa"
  nofw:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --nofw
    doc: |
      --norc      do not align to forward/reverse-complement reference strand
  phred33-quals:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --phred33-quals
    doc: "input quals are Phred+33 (default)"
  k:
    type: int?
    inputBinding:
      position: 1
      prefix: -k
    doc: |
      <int>           report up to <int> good alignments per read (default: 1)
  filelist_mates:
    type: File[]?
    inputBinding:
      itemSeparator: ','
      position: 10

  m:
    type: int?
    inputBinding:
      position: 1
      prefix: -m
    doc: |
      <int>           suppress all alignments if > <int> exist (def: no limit)
  l:
    type: int?
    inputBinding:
      position: 1
      prefix: -l
    doc: |
      --seedlen <int> seed length for -n (default: 28)
  filelist:
    type:
      type: array
      items: File
    inputBinding:
      itemSeparator: ','
      position: 9

    doc: |
      {-1 <m1> -2 <m2> | --12 <r> | <s>}
      <m1>    Comma-separated list of files containing upstream mates (or the
            sequences themselves, if -c is set) paired with mates in <m2>
      <m2>    Comma-separated list of files containing downstream mates (or the
            sequences themselves if -c is set) paired with mates in <m1>
      <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
            a mixture of paired and unpaired.  Specify "-"for stdin.
      <s>     Comma-separated list of files containing unpaired reads, or the
            sequences themselves, if -c is set.  Specify "-"for stdin.
  n:
    type: int?
    inputBinding:
      position: 1
      prefix: -n
    doc: |
      --seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
  q:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -q
    doc: "query input files are FASTQ .fq/.fastq (default)"
  p:
    type: int?
    inputBinding:
      position: 1
      prefix: -p
    doc: |
      --threads <int> number of alignment threads to launch (default: 1)
  s:
    type: int?
    inputBinding:
      position: 1
      prefix: -s
    doc: |
      --skip <int>    skip the first <int> reads/pairs in the input
  r:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -r
    doc: "query input files are raw one-sequence-per-line"
  u:
    type: int?
    inputBinding:
      position: 1
      prefix: -u
    doc: |
      --qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
  t:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -t
    doc: |
      --time          print wall-clock time taken by search phases
  refidx:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --refidx
    doc: "refer to ref. seqs by 0-based index rather than name"
  v:
    type: int?
    inputBinding:
      position: 1
      prefix: -v
    doc: "<int>           report end-to-end hits w/ <=v mismatches; ignore qualities"
  integer-quals:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --integer-quals
    doc: "qualities are given as space-separated integers (not ASCII)"
  solexa-quals:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --solexa-quals
    doc: "input quals are from GA Pipeline ver. < 1.3"
  snpfrac:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --snpfrac
    doc: |
      <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
  refout:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --refout
    doc: "write alignments to files refXXXXX.map, 1 map per reference"
  nomaqround:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --nomaqround
    doc: "disable Maq-like quality rounding for -n (nearest 10 <= 30)"
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.filename)

  output_bowtie_log:
    type: File
    outputBinding:
      glob: $(inputs.filename + '.log')

baseCommand:
- bowtie

arguments:
- valueFrom: $('2> ' + inputs.filename + '.log')
  position: 100000
  shellQuote: false

$namespaces:
  schema: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

schema:mainEntity:
#  $import: https://scidap.com/description/tools/bowtie.yaml
  class: schema:SoftwareSourceCode
  schema:name: bowtie
  schema:about: 'Bowtie is an ultrafast, memory-efficient short read aligner. It aligns
    short DNA sequences (reads) to the human genome at a rate of over 25 million 35-bp
    reads per hour. Bowtie indexes the genome with a Burrows-Wheeler index to keep
    its memory footprint small: typically about 2.2 GB for the human genome (2.9 GB
    for paired-end).

    '
  schema:url: http://bowtie-bio.sourceforge.net
  schema:codeRepository: https://github.com/BenLangmead/bowtie.git

  schema:license:
  - https://opensource.org/licenses/GPL-3.0

  schema:targetProduct:
    class: schema:SoftwareApplication
    schema:softwareVersion: 1.1.2
    schema:applicationCategory: commandline tool
  schema:programmingLanguage: C++
  schema:publication:
  - class: schema:ScholarlyArticle
    id: https://doi.org/10.1186/gb-2009-10-3-r25

schema:author:
  class: schema:Person
  schema:name: Andrey Kartashov
  schema:email: mailto:Andrey.Kartashov@cchmc.org
  schema:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  schema:worksFor:
  - class: schema:Organization
    schema:name: Cincinnati Children's Hospital Medical Center
    schema:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    schema:department:
    - class: schema:Organization
      schema:name: Barski Lab
doc: "bowtie.cwl is developed for CWL consortium"

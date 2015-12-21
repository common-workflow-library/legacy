#!/usr/bin/env cwl-runner

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  adms: http://www.w3.org/ns/adms#
  dcat: http://www.w3.org/ns/dcat#

$schemas:
- http://schema.rdfs.org/all.rdf
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

adms:includedAsset:
  doap:name: "bowtie"
  doap:description: ""
  doap:homepage: "http://bowtie-bio.sourceforge.net"
  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/BenLangmead/bowtie.git"
  doap:release:
  - class: doap:Version
    doap:revision: "1.1.2"
  doap:license: "GPL"
  doap:category: "commandline tool"
  doap:programming-language: "C++"
  foaf:publications:
  - id: urn:pmid:0
    foaf:title: ""
    foaf:homepage: ""
  doap:developer:
  - class: foaf:Organization
    foaf:name: ""
    foaf:member:
    - class: foaf:Person
      foaf:name: ""
      foaf:mbox: ""
      foaf:fundedBy: ""

label: "bowtie.cwl"
doap:name: "bowtie.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/bowtie.cwl"

dct:isPartOf:
  doap:name: "CWL Workflows"
  doap:homepage: "http://commonwl.org/"
  doap:license: "Apache2"

  doap:implements:
  - class: doap:Specification
    doap:homepage: "http://common-workflow-language.github.io/draft-3/"

  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/common-workflow-language/workflows"

  dct:creator:
  - class: foaf:Organization
    foaf:name: "Curoverse"
    foaf:member:
    - class: foaf:Person
      id: "http://orcid.org/0000-0003-3566-7705"
      foaf:name: "Peter Amstutz"
      foaf:mbox: "mailto:peter.amstutz@curoverse.com"
  - class: foaf:Organization
    foaf:name: "Seven Bridges Genomics"
    foaf:member:
    - class: foaf:Person
      id: "mailto:nebojsa.tijanic@sbgenomics.com"
      foaf:name: "Nebojša Tijanić"
      foaf:mbox: "mailto:nebojsa.tijanic@sbgenomics.com"

  dct:contributor:
  - class: foaf:Organization
    foaf:name: "Seven Bridges Genomics"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Luka Stojanovic"
      foaf:mbox: "mailto:luka.stojanovic@sbgenomics.com"
  - class: foaf:Organization
    foaf:name: "Galaxy Project, Pennsylvania State University"
    foaf:member:
    - class: foaf:Person
      foaf:name: "John Chilton"
      foaf:mbox: "mailto:jmchilton@gmail.com"
  - class: foaf:Organization
    foaf:name: "University of California, Davis"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Michael R. Crusoe"
      foaf:mbox: "mailto:crusoe@ucdavis.edu"
  - class: foaf:Organization
    foaf:name: "Institut Pasteur"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Hervé Ménager"
      foaf:mbox: "mailto:herve.menager@gmail.com"
  - class: foaf:Organization
    foaf:name: "BioDatomics"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Maxim Mikheev"
      foaf:mbox: "mailto:mikhmv@biodatomics.com"
  - class: foaf:Organization
    foaf:name: "University of Manchester"
    foaf:member:
    - class: foaf:Person
      foaf:name: "Stian Soiland-Reyes"
      foaf:mbox: "mailto:soiland-reyes@cs.manchester.ac.uk"

doap:maintainer:
- class: foaf:Organization
  foaf:name: "Barski Lab, Cincinnati Children's Hospital Medical Center"
  foaf:member:
  - class: foaf:Person
    id: "http://orcid.org/0000-0001-9102-5681"
    foaf:openid: "http://orcid.org/0000-0001-9102-5681"
    foaf:name: "Andrey Kartashov"
    foaf:mbox: "mailto:Andrey.Kartashov@cchmc.org"


requirements:
  - $import: envvar-global.cwl
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    #dockerImageId: scidap/bowtie:v1.1.2 #not yet ready
    dockerPull: scidap/bowtie:v1.1.2
    dockerFile: |
      #################################################################
      # Dockerfile
      #
      # Software:         bowtie
      # Software Version: 1.1.2
      # Description:      Bowtie image for SciDAP
      # Website:          http://bowtie-bio.sourceforge.net, http://scidap.com/
      # Provides:         bowtie
      # Base Image:       scidap/scidap:v0.0.1
      # Build Cmd:        docker build --rm -t scidap/bowtie:v1.1.2 .
      # Pull Cmd:         docker pull scidap/bowtie:v1.1.2
      # Run Cmd:          docker run --rm scidap/bowtie:v1.1.2 bowtie
      #################################################################

      ### Base Image
      FROM scidap/scidap:v0.0.1
      MAINTAINER Andrey V Kartashov "porter@porter.st"
      ENV DEBIAN_FRONTEND noninteractive

      ################## BEGIN INSTALLATION ######################

      WORKDIR /tmp

      ### Installing bowtie

      ENV VERSION 1.1.2
      ENV NAME bowtie
      ENV URL "https://github.com/BenLangmead/bowtie/archive/v${VERSION}.tar.gz"

      RUN wget -q -O - $URL | tar -zxv && \
          cd ${NAME}-${VERSION} && \
          make -j 4 && \
          cd .. && \
          cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
          cp ./${NAME}-${VERSION}/${NAME}-* /usr/local/bin/ && \
          strip /usr/local/bin/*; true && \
          rm -rf ./${NAME}-${VERSION}/


inputs:
  - id: '#stdoutfile'
    type: string

  - id: '#q'
    type:
      - 'null'
      - boolean
    description: "query input files are FASTQ .fq/.fastq (default)\n"
    inputBinding:
      position: 1
      prefix: '-q'

  - id: '#f'
    type:
      - 'null'
      - boolean
    description: "query input files are (multi-)FASTA .fa/.mfa\n"
    inputBinding:
      position: 1
      prefix: '-f'

  - id: '#r'
    type:
      - 'null'
      - boolean
    description: "query input files are raw one-sequence-per-line\n"
    inputBinding:
      position: 1
      prefix: '-r'

  - id: '#c'
    type:
      - 'null'
      - boolean
    description: "query sequences given on cmd line (as <mates>, <singles>)\n"
    inputBinding:
      position: 1
      prefix: '-c'

  - id: '#C'
    type:
      - 'null'
      - boolean
    description: "reads and index are in colorspace\n"
    inputBinding:
      position: 1
      prefix: '-C'
  - id: '#Q'
    type:
      - 'null'
      - File
    description: >
      --quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
    inputBinding:
      position: 1
      prefix: '-Q'
  - id: '#Q1'
    type:
      - 'null'
      - boolean
    description: |
      --Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
    inputBinding:
      position: 1
      prefix: '--Q1'
  - id: '#s'
    type:
      - 'null'
      - int
    description: |
      --skip <int>    skip the first <int> reads/pairs in the input
    inputBinding:
      position: 1
      prefix: '-s'
  - id: '#u'
    type:
      - 'null'
      - int
    description: |
      --qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
    inputBinding:
      position: 1
      prefix: '-u'
  - id: '#5'
    type:
      - 'null'
      - int
    description: |
      --trim5 <int>   trim <int> bases from 5' (left) end of reads
    inputBinding:
      position: 1
      prefix: '-5'
  - id: '#3'
    type:
      - 'null'
      - int
    description: |
      --trim3 <int>   trim <int> bases from 3' (right) end of reads
    inputBinding:
      position: 1
      prefix: '-3'
  - id: '#phred33-quals'
    type:
      - 'null'
      - boolean
    description: "input quals are Phred+33 (default)\n"
    inputBinding:
      position: 1
      prefix: '--phred33-quals'
  - id: '#phred64-quals'
    type:
      - 'null'
      - boolean
    description: "input quals are Phred+64 (same as --solexa1.3-quals)\n"
    inputBinding:
      position: 1
      prefix: '--phred64-quals'
  - id: '#solexa-quals'
    type:
      - 'null'
      - boolean
    description: "input quals are from GA Pipeline ver. < 1.3\n"
    inputBinding:
      position: 1
      prefix: '--solexa-quals'
  - id: '#solexa1.3-quals'
    type:
      - 'null'
      - boolean
    description: "input quals are from GA Pipeline ver. >= 1.3\n"
    inputBinding:
      position: 1
      prefix: '--solexa1.3-quals'
  - id: '#integer-quals'
    type:
      - 'null'
      - boolean
    description: "qualities are given as space-separated integers (not ASCII)\n"
    inputBinding:
      position: 1
      prefix: '--integer-quals'
  - id: '#large-index'
    type:
      - 'null'
      - boolean
    description: "force usage of a 'large' index, even if a small one is present\nAlignment:\n"
    inputBinding:
      position: 1
      prefix: '--large-index'
  - id: '#v'
    type:
      - 'null'
      - int
    description: >
      <int>           report end-to-end hits w/ <=v mismatches; ignore qualities

      or
    inputBinding:
      position: 1
      prefix: '-v'
  - id: '#n'
    type:
      - 'null'
      - int
    description: |
      --seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
    inputBinding:
      position: 1
      prefix: '-n'
  - id: '#e'
    type:
      - 'null'
      - int
    description: >
      --maqerr <int>  max sum of mismatch quals across alignment for -n (def:
      70)
    inputBinding:
      position: 1
      prefix: '-e'
  - id: '#l'
    type:
      - 'null'
      - int
    description: |
      --seedlen <int> seed length for -n (default: 28)
    inputBinding:
      position: 1
      prefix: '-l'
  - id: '#nomaqround'
    type:
      - 'null'
      - boolean
    description: "disable Maq-like quality rounding for -n (nearest 10 <= 30)\n"
    inputBinding:
      position: 1
      prefix: '--nomaqround'
  - id: '#I'
    type:
      - 'null'
      - int
    description: |
      --minins <int>  minimum insert size for paired-end alignment (default: 0)
    inputBinding:
      position: 1
      prefix: '-I'
  - id: '#X'
    type:
      - 'null'
      - int
    description: >
      --maxins <int>  maximum insert size for paired-end alignment (default:
      250)
    inputBinding:
      position: 1
      prefix: '-X'
  - id: '#fr'
    type:
      - 'null'
      - boolean
    description: |
      --rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
    inputBinding:
      position: 1
      prefix: '--fr'
  - id: '#nofw'
    type:
      - 'null'
      - boolean
    description: |
      --norc      do not align to forward/reverse-complement reference strand
    inputBinding:
      position: 1
      prefix: '--nofw'
  - id: '#maxbts'
    type:
      - 'null'
      - int
    description: |
      <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
    inputBinding:
      position: 1
      prefix: '--maxbts'
  - id: '#pairtries'
    type:
      - 'null'
      - int
    description: |
      <int>  max # attempts to find mate for anchor hit (default: 100)
    inputBinding:
      position: 1
      prefix: '--pairtries'
  - id: '#y'
    type:
      - 'null'
      - boolean
    description: >
      --tryhard       try hard to find valid alignments, at the expense of speed
    inputBinding:
      position: 1
      prefix: '-y'
  - id: '#chunkmbs'
    type:
      - 'null'
      - int
    description: |
      <int>   max megabytes of RAM for best-first search frames (def: 64)
      Reporting:
    inputBinding:
      position: 1
      prefix: '--chunkmbs'
  - id: '#k'
    type:
      - 'null'
      - int
    description: |
      <int>           report up to <int> good alignments per read (default: 1)
    inputBinding:
      position: 1
      prefix: '-k'
  - id: '#a'
    type:
      - 'null'
      - boolean
    description: |
      --all           report all alignments per read (much slower than low -k)
    inputBinding:
      position: 1
      prefix: '-a'
  - id: '#m'
    type:
      - 'null'
      - int
    description: |
      <int>           suppress all alignments if > <int> exist (def: no limit)
    inputBinding:
      position: 1
      prefix: '-m'
  - id: '#M'
    type:
      - 'null'
      - int
    description: >
      <int>           like -m, but reports 1 random hit (MAPQ=0); requires
      --best
    inputBinding:
      position: 1
      prefix: '-M'
  - id: '#best'
    type:
      - 'null'
      - boolean
    description: "hits guaranteed best stratum; ties broken by quality\n"
    inputBinding:
      position: 1
      prefix: '--best'
  - id: '#strata'
    type:
      - 'null'
      - boolean
    description: "hits in sub-optimal strata aren't reported (requires --best)\nOutput:\n"
    inputBinding:
      position: 1
      prefix: '--strata'
  - id: '#t'
    type:
      - 'null'
      - boolean
    description: |
      --time          print wall-clock time taken by search phases
    inputBinding:
      position: 1
      prefix: '-t'
  - id: '#B'
    type:
      - 'null'
      - int
    description: |
      --offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
    inputBinding:
      position: 1
      prefix: '-B'
  - id: '#quiet'
    type:
      - 'null'
      - boolean
    description: "print nothing but the alignments\n"
    inputBinding:
      position: 1
      prefix: '--quiet'
  - id: '#refout'
    type:
      - 'null'
      - boolean
    description: "write alignments to files refXXXXX.map, 1 map per reference\n"
    inputBinding:
      position: 1
      prefix: '--refout'
  - id: '#refidx'
    type:
      - 'null'
      - boolean
    description: "refer to ref. seqs by 0-based index rather than name\n"
    inputBinding:
      position: 1
      prefix: '--refidx'
  - id: '#al'
    type:
      - 'null'
      - boolean
    description: |
      <fname>       write aligned reads/pairs to file(s) <fname>
    inputBinding:
      position: 1
      prefix: '--al'
  - id: '#un'
    type:
      - 'null'
      - boolean
    description: |
      <fname>       write unaligned reads/pairs to file(s) <fname>
    inputBinding:
      position: 1
      prefix: '--un'
  - id: '#max'
    type:
      - 'null'
      - boolean
    description: |
      <fname>      write reads/pairs over -m limit to file(s) <fname>
    inputBinding:
      position: 1
      prefix: '--max'
  - id: '#suppress'
    type:
      - 'null'
      - boolean
    description: |
      <cols>  suppresses given columns (comma-delim'ed) in default output
    inputBinding:
      position: 1
      prefix: '--suppress'
  - id: '#fullref'
    type:
      - 'null'
      - boolean
    description: "write entire ref name (default: only up to 1st space)\nColorspace:\n"
    inputBinding:
      position: 1
      prefix: '--fullref'
  - id: '#snpphred'
    type:
      - 'null'
      - int
    description: |
      <int>   Phred penalty for SNP when decoding colorspace (def: 30)
      or
    inputBinding:
      position: 1
      prefix: '--snpphred'
  - id: '#snpfrac'
    type:
      - 'null'
      - boolean
    description: |
      <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
    inputBinding:
      position: 1
      prefix: '--snpfrac'
  - id: '#col-cseq'
    type:
      - 'null'
      - boolean
    description: "print aligned colorspace seqs as colors, not decoded bases\n"
    inputBinding:
      position: 1
      prefix: '--col-cseq'
  - id: '#col-cqual'
    type:
      - 'null'
      - boolean
    description: "print original colorspace quals, not decoded quals\n"
    inputBinding:
      position: 1
      prefix: '--col-cqual'
  - id: '#col-keepends'
    type:
      - 'null'
      - boolean
    description: "keep nucleotides at extreme ends of decoded alignment\nSAM:\n"
    inputBinding:
      position: 1
      prefix: '--col-keepends'
  - id: '#S'
    type:
      - 'null'
      - boolean
    description: |
      --sam           write hits in SAM format
    inputBinding:
      position: 1
      prefix: '-S'
  - id: '#mapq'
    type:
      - 'null'
      - int
    description: |
      <int>       default mapping quality (MAPQ) to print for SAM alignments
    inputBinding:
      position: 1
      prefix: '--mapq'
  - id: '#sam-nohead'
    type:
      - 'null'
      - boolean
    description: "supppress header lines (starting with @) for SAM output\n"
    inputBinding:
      position: 1
      prefix: '--sam-nohead'
  - id: '#sam-nosq'
    type:
      - 'null'
      - boolean
    description: "supppress @SQ header lines for SAM output\n"
    inputBinding:
      position: 1
      prefix: '--sam-nosq'
  - id: '#sam-RG'
    type:
      - 'null'
      - string
    description: |
      <text>    add <text> (usually "lab=value") to @RG line of SAM header
      Performance:
    inputBinding:
      position: 1
      prefix: '--sam-RG'
  - id: '#o'
    type:
      - 'null'
      - int
    description: |
      --offrate <int> override offrate of index; must be >= index's offrate
    inputBinding:
      position: 1
      prefix: '-o'
  - id: '#p'
    type:
      - 'null'
      - int
    description: |
      --threads <int> number of alignment threads to launch (default: 1)
    inputBinding:
      position: 1
      prefix: '-p'
  - id: '#mm'
    type:
      - 'null'
      - boolean
    description: "use memory-mapped I/O for index; many 'bowtie's can share\n"
    inputBinding:
      position: 1
      prefix: '--mm'
  - id: '#shmem'
    type:
      - 'null'
      - boolean
    description: "use shared mem for index; many 'bowtie's can share\nOther:\n"
    inputBinding:
      position: 1
      prefix: '--shmem'
  - id: '#seed'
    type:
      - 'null'
      - int
    description: |
      <int>       seed for random number generator
    inputBinding:
      position: 1
      prefix: '--seed'
  - id: '#verbose'
    type:
      - 'null'
      - boolean
    description: "verbose output (for debugging)\n"
    inputBinding:
      position: 1
      prefix: '--verbose'

outputs:
  - id: '#stdoutfile'
    type: File
    outputBinding:
      glob: $(inputs.stdoutfile)

stdout: $(inputs.stdoutfile)

baseCommand:
  - bowtie

description: |
  bowtie.cwl is developed for CWL consortium

  Usage: 
  bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]

    <m1>    Comma-separated list of files containing upstream mates (or the
            sequences themselves, if -c is set) paired with mates in <m2>
    <m2>    Comma-separated list of files containing downstream mates (or the
            sequences themselves if -c is set) paired with mates in <m1>
    <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
            a mixture of paired and unpaired.  Specify "-"for stdin.
    <s>     Comma-separated list of files containing unpaired reads, or the
            sequences themselves, if -c is set.  Specify "-"for stdin.
    <hit>   File to write hits to (default: stdout)
  Input:
    -q                 query input files are FASTQ .fq/.fastq (default)
    -f                 query input files are (multi-)FASTA .fa/.mfa
    -r                 query input files are raw one-sequence-per-line
    -c                 query sequences given on cmd line (as <mates>, <singles>)
    -C                 reads and index are in colorspace
    -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
    --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
    -s/--skip <int>    skip the first <int> reads/pairs in the input
    -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
    -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
    -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads
    --phred33-quals    input quals are Phred+33 (default)
    --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)
    --solexa-quals     input quals are from GA Pipeline ver. < 1.3
    --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
    --integer-quals    qualities are given as space-separated integers (not ASCII)
    --large-index      force usage of a 'large' index, even if a small one is present
  Alignment:
    -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
      or
    -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
    -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
    -l/--seedlen <int> seed length for -n (default: 28)
    --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
    -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)
    -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)
    --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
    --nofw/--norc      do not align to forward/reverse-complement reference strand
    --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
    --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)
    -y/--tryhard       try hard to find valid alignments, at the expense of speed
    --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
  Reporting:
    -k <int>           report up to <int> good alignments per read (default: 1)
    -a/--all           report all alignments per read (much slower than low -k)
    -m <int>           suppress all alignments if > <int> exist (def: no limit)
    -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
    --best             hits guaranteed best stratum; ties broken by quality
    --strata           hits in sub-optimal strata aren't reported (requires --best)
  Output:
    -t/--time          print wall-clock time taken by search phases
    -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
    --quiet            print nothing but the alignments
    --refout           write alignments to files refXXXXX.map, 1 map per reference
    --refidx           refer to ref. seqs by 0-based index rather than name
    --al <fname>       write aligned reads/pairs to file(s) <fname>
    --un <fname>       write unaligned reads/pairs to file(s) <fname>
    --max <fname>      write reads/pairs over -m limit to file(s) <fname>
    --suppress <cols>  suppresses given columns (comma-delim'ed) in default output
    --fullref          write entire ref name (default: only up to 1st space)
  Colorspace:
    --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)
       or
    --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
    --col-cseq         print aligned colorspace seqs as colors, not decoded bases
    --col-cqual        print original colorspace quals, not decoded quals
    --col-keepends     keep nucleotides at extreme ends of decoded alignment
  SAM:
    -S/--sam           write hits in SAM format
    --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments
    --sam-nohead       supppress header lines (starting with @) for SAM output
    --sam-nosq         supppress @SQ header lines for SAM output
    --sam-RG <text>    add <text> (usually "lab=value") to @RG line of SAM header
  Performance:
    -o/--offrate <int> override offrate of index; must be >= index's offrate
    -p/--threads <int> number of alignment threads to launch (default: 1)
    --mm               use memory-mapped I/O for index; many 'bowtie's can share
    --shmem            use shared mem for index; many 'bowtie's can share
  Other:
    --seed <int>       seed for random number generator
    --verbose          verbose output (for debugging)
    --version          print version information and quit
    -h/--help          print this usage message


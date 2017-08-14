#!/usr/bin/env cwl-runner

# Mantainer: alejandro.barrera@duke.edu
# Partially Auto generated with clihp (https://github.com/portah/clihp, developed by Andrey.Kartashov@cchmc.org)
# Developed for GGR project (https://github.com/Duke-GCB/GGR-cwl)

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: deeptools-docker.yml
- class: InlineJavascriptRequirement

inputs:
  verbose:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --verbose
    doc: |
      --verbose
      Set to see processing messages. (default: False)
  binSize:
    type: int?
    inputBinding:
      position: 1
      prefix: --binSize
    doc: |
      INT bp
      Size of the bins, in bases, for the output of the
      bigwig/bedgraph file. (default: 50)
  MNase:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --MNase
    doc: |
      Determine nucleosome positions from MNase-seq data.
      Only 3 nucleotides at the center of each fragment are
      counted. The fragment ends are defined by the two mate
      reads. Only fragment lengthsbetween 130 - 200 bp are
      considered to avoid dinucleosomes or other
      artifacts.*NOTE*: Requires paired-end data. A bin size
      of 1 is recommended. (default: False)
  ignoreDuplicates:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --ignoreDuplicates
    doc: |
      If set, reads that have the same orientation and start
      position will be considered only once. If reads are
      paired, the mate's position also has to coincide to
      ignore a read. (default: False)
  numberOfProcessors:
    type: int?
    inputBinding:
      position: 1
      prefix: --numberOfProcessors
    doc: |
      INT
      Number of processors to use. Type "max/2" to use half
      the maximum number of processors or "max" to use all
      available processors. (default: max/2)
  ignoreForNormalization:
    type: string?
    inputBinding:
      position: 1
      prefix: --ignoreForNormalization
    doc: |
      --ignoreForNormalization chrX chrM. (default: None)
      A list of space-delimited chromosome names containing
      those chromosomes that should be excluded for
      computing the normalization. This is useful when
      considering samples with unequal coverage across
      chromosomes, like male samples. An usage examples is
  outFileName:
#-------------------------------------
#--- Output formatting arguments -----
#-------------------------------------
    type: string?
    doc: |
      FILENAME
      Output file name. (default: input BAM filename with bigwig [*.bw] or bedgraph [*.bdg] extension.)
  smoothLength:
    type: int?
    inputBinding:
      position: 1
      prefix: --smoothLength
    doc: |
      INT bp
      The smooth length defines a window, larger than the
      binSize, to average the number of reads. For example,
      if the --binSize is set to 20 and the --smoothLength
      is set to 60, then, for each bin, the average of the
      bin and its left and right neighbors is considered.
      Any value smaller than --binSize will be ignored and
      no smoothing will be applied. (default: None)
      Read processing options:
  version:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --version
    doc: show program's version number and exit
  extendReads:
#-----------------------------------------
#------- Read processing options ---------
#-----------------------------------------
    type: int?
    inputBinding:
      position: 1
      prefix: --extendReads
    doc: |
      INT bp
      This parameter allows the extension of reads to
      fragment size. If set, each read is extended, without
      exception. *NOTE*: This feature is generally NOT
      recommended for spliced-read data, such as RNA-seq, as
      it would extend reads over skipped regions. *Single-
      end*: Requires a user specified value for the final
      fragment length. Reads that already exceed this
      fragment length will not be extended. *Paired-end*:
      Reads with mates are always extended to match the
      fragment size defined by the two read mates. Unmated
      reads, mate reads that map too far apart (>4x fragment
      length) or even map to different chromosomes are
      treated like single-end reads. The input of a fragment
      length value is optional. If no value is specified, it
      is estimated from the data (mean of the fragment size
      of all mate reads). (default: False)
  centerReads:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --centerReads
    doc: |
      By adding this option, reads are centered with respect
      to the fragment length. For paired-end data, the read
      is centered at the fragment length defined by the two
      ends of the fragment. For single-end data, the given
      fragment length is used. This option is useful to get
      a sharper signal around enriched regions. (default:
      False)
  samFlagExclude:
    type: int?
    inputBinding:
      position: 1
      prefix: --samFlagExclude
    doc: |
      INT
      Exclude reads based on the SAM flag. For example, to
      get only reads that map to the forward strand, use
      --samFlagExclude 16, where 16 is the SAM flag for
      reads that map to the reverse strand. (default: None)
  samFlagInclude:
    type: int?
    inputBinding:
      position: 1
      prefix: --samFlagInclude
    doc: |
      INT
      Include reads based on the SAM flag. For example, to
      get only reads that are the first mate, use a flag of
      64. This is useful to count properly paired reads only
      once, as otherwise the second mate will be also
      considered for the coverage. (default: None)
  filterRNAstrand:
    type: string?
    inputBinding:
      position: 1
      prefix: --filterRNAstrand
    doc: |
      {forward,reverse}
      Selects RNA-seq reads (single-end or paired-end) in
      the given strand. (default: None)
  scaleFactor:
#---------------------------------
#-----  Processing arguments -----
#---------------------------------
    type: float?
    inputBinding:
      position: 1
      prefix: --scaleFactor
    doc: |
      SCALEFACTOR
      Indicate a number that you would like to use. When
      used in combination with --normalizeTo1x or
      --normalizeUsingRPKM, the computed scaling factor will
      be multiplied by the given scale factor. (default:
      1.0)
  skipNonCoveredRegions:
    type: string?
    inputBinding:
      position: 1
      prefix: --skipNonCoveredRegions
    doc: |
      --skipNonCoveredRegions, --skipNAs
      This parameter determines if non-covered regions
      (regions without overlapping reads) in a BAM file
      should be skipped. The default is to treat those
      regions as having a value of zero. The decision to
      skip non-covered regions depends on the interpretation
      of the data. Non-covered regions may represent, for
      example, repetitive regions that should be skipped.
      (default: False)
  outFileFormat:
    type: string
    default: bigwig
    inputBinding:
      position: 1
      prefix: --outFileFormat
    doc: |
      {bigwig,bedgraph}, -of {bigwig,bedgraph}
      Output file type. Either "bigwig" or "bedgraph".
      (default: bigwig)
  output_suffix:
    type: string?
    doc: Suffix used for output file (input BAM filename + suffix)
  region:
    type: string?
    inputBinding:
      position: 1
      prefix: --region
    doc: |
      CHR:START:END
      Region of the genome to limit the operation to - this
      is useful when testing parameters to reduce the
      computing time. The format is chr:start:end, for
      example --region chr10 or --region
      chr10:456700:891000. (default: None)
  normalizeUsingRPKM:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --normalizeUsingRPKM
    doc: |
      Use Reads Per Kilobase per Million reads to normalize
      the number of reads per bin. The formula is: RPKM (per
      bin) = number of reads per bin / ( number of mapped
      reads (in millions) * bin length (kb) ). Each read is
      considered independently,if you want to only count
      either of the mate pairs inpaired-end data, use the
      --samFlag option. (default: False)
  normalizeTo1x:
#-----------------------------------------
#-- Read coverage normalization options --
#-----------------------------------------
    type: string?
    inputBinding:
      position: 1
      prefix: --normalizeTo1x
    doc: |
      EFFECTIVE GENOME SIZE LENGTH
      Report read coverage normalized to 1x sequencing depth
      (also known as Reads Per Genomic Content (RPGC)).
      Sequencing depth is defined as: (total number of
      mapped reads * fragment length) / effective genome
      size. The scaling factor used is the inverse of the
      sequencing depth computed for the sample to match the
      1x coverage. To use this option, the effective genome
      size has to be indicated after the option. The
      effective genome size is the portion of the genome
      that is mappable. Large fractions of the genome are
      stretches of NNNN that should be discarded. Also, if
      repetitive regions were not included in the mapping of
      reads, the effective genome size needs to be adjusted
      accordingly. Common values are: mm9: 2,150,570,000;
      hg19:2,451,960,000; dm3:121,400,000 and
      ce10:93,260,000. See Table 2 of http://www.plosone.org
      /article/info:doi/10.1371/journal.pone.0030377 or http
      ://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_
      T1.html for several effective genome sizes. (default:
      None)
  blackListFileName:
    type: File?
    inputBinding:
      position: 1
      prefix: --blackListFileName
    doc: |
      BED file
      A BED file containing regions that should be excluded
      from all analyses. Currently this works by rejecting
      genomic chunks that happen to overlap an entry.
      Consequently, for BAM files, if a read partially
      overlaps a blacklisted region or a fragment spans over
      it, then the read/fragment might still be considered.
      (default: None)
  bam:
    type: File
    secondaryFiles: $(self.path + '.bai')
    inputBinding:
      position: 1
      prefix: --bam
    doc: 'BAM file to process '
  minMappingQuality:
    type: int?
    inputBinding:
      position: 1
      prefix: --minMappingQuality
    doc: |
      INT
      If set, only reads that have a mapping quality score
      of at least this are considered. (default: None)
outputs:
  output_bam_coverage:
    type: File
    outputBinding:
      glob: ${ if (inputs.outFileName) return inputs.outFileName; if (inputs.output_suffix)
        return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") +
        inputs.output_suffix; if (inputs.outFileFormat == "bedgraph") return inputs.bam.path.replace(/^.*[\\\/]/,
        "").replace(/\.[^/.]+$/, "") + ".bdg"; return inputs.bam.path.replace(/^.*[\\\/]/,
        "").replace(/\.[^/.]+$/, "") + ".bw"; }
baseCommand: bamCoverage
arguments:
- valueFrom: ${ if (inputs.outFileName) return inputs.outFileName; if (inputs.output_suffix)
    return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + inputs.output_suffix;
    if (inputs.outFileFormat == "bedgraph") return inputs.bam.path.replace(/^.*[\\\/]/,
    "").replace(/\.[^/.]+$/, "") + ".bdg"; return inputs.bam.path.replace(/^.*[\\\/]/,
    "").replace(/\.[^/.]+$/, "") + ".bw"; }
  prefix: --outFileName
  position: 3
doc: |
  usage: An example usage is:$ bamCoverage -b reads.bam -o coverage.bw


  This tool takes an alignment of reads or fragments as input (BAM file) and
  generates a coverage track (bigWig or bedGraph) as output. The coverage is
  calculated as the number of reads per bin, where bins are short consecutive
  counting windows of a defined size. It is possible to extended the length of
  the reads to better reflect the actual fragment length. *bamCoverage* offers
  normalization by scaling factor, Reads Per Kilobase per Million mapped reads
  (RPKM), and 1x depth (reads per genome coverage, RPGC).

  Required arguments:
    --bam BAM file, -b BAM file
                          BAM file to process (default: None)

  Output:
    --outFileName FILENAME, -o FILENAME
                          Output file name. (default: None)
    --outFileFormat {bigwig,bedgraph}, -of {bigwig,bedgraph}
                          Output file type. Either "bigwig" or "bedgraph".
                          (default: bigwig)

  Optional arguments:
    --help, -h            show this help message and exit
    --scaleFactor SCALEFACTOR
                          Indicate a number that you would like to use. When
                          used in combination with --normalizeTo1x or
                          --normalizeUsingRPKM, the computed scaling factor will
                          be multiplied by the given scale factor. (default:
                          1.0)
    --MNase               Determine nucleosome positions from MNase-seq data.
                          Only 3 nucleotides at the center of each fragment are
                          counted. The fragment ends are defined by the two mate
                          reads. Only fragment lengthsbetween 130 - 200 bp are
                          considered to avoid dinucleosomes or other
                          artifacts.*NOTE*: Requires paired-end data. A bin size
                          of 1 is recommended. (default: False)
    --filterRNAstrand {forward,reverse}
                          Selects RNA-seq reads (single-end or paired-end) in
                          the given strand. (default: None)
    --version             show program's version number and exit
    --binSize INT bp, -bs INT bp
                          Size of the bins, in bases, for the output of the
                          bigwig/bedgraph file. (default: 50)
    --region CHR:START:END, -r CHR:START:END
                          Region of the genome to limit the operation to - this
                          is useful when testing parameters to reduce the
                          computing time. The format is chr:start:end, for
                          example --region chr10 or --region
                          chr10:456700:891000. (default: None)
    --blackListFileName BED file, -bl BED file
                          A BED file containing regions that should be excluded
                          from all analyses. Currently this works by rejecting
                          genomic chunks that happen to overlap an entry.
                          Consequently, for BAM files, if a read partially
                          overlaps a blacklisted region or a fragment spans over
                          it, then the read/fragment might still be considered.
                          (default: None)
    --numberOfProcessors INT, -p INT
                          Number of processors to use. Type "max/2" to use half
                          the maximum number of processors or "max" to use all
                          available processors. (default: max/2)
    --verbose, -v         Set to see processing messages. (default: False)

  Read coverage normalization options:
    --normalizeTo1x EFFECTIVE GENOME SIZE LENGTH
                          Report read coverage normalized to 1x sequencing depth
                          (also known as Reads Per Genomic Content (RPGC)).
                          Sequencing depth is defined as: (total number of
                          mapped reads * fragment length) / effective genome
                          size. The scaling factor used is the inverse of the
                          sequencing depth computed for the sample to match the
                          1x coverage. To use this option, the effective genome
                          size has to be indicated after the option. The
                          effective genome size is the portion of the genome
                          that is mappable. Large fractions of the genome are
                          stretches of NNNN that should be discarded. Also, if
                          repetitive regions were not included in the mapping of
                          reads, the effective genome size needs to be adjusted
                          accordingly. Common values are: mm9: 2,150,570,000;
                          hg19:2,451,960,000; dm3:121,400,000 and
                          ce10:93,260,000. See Table 2 of http://www.plosone.org
                          /article/info:doi/10.1371/journal.pone.0030377 or http
                          ://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_
                          T1.html for several effective genome sizes. (default:
                          None)
    --normalizeUsingRPKM  Use Reads Per Kilobase per Million reads to normalize
                          the number of reads per bin. The formula is: RPKM (per
                          bin) = number of reads per bin / ( number of mapped
                          reads (in millions) * bin length (kb) ). Each read is
                          considered independently,if you want to only count
                          either of the mate pairs inpaired-end data, use the
                          --samFlag option. (default: False)
    --ignoreForNormalization IGNOREFORNORMALIZATION [IGNOREFORNORMALIZATION ...],
  -ignore IGNOREFORNORMALIZATION [IGNOREFORNORMALIZATION ...]
                          A list of space-delimited chromosome names containing
                          those chromosomes that should be excluded for
                          computing the normalization. This is useful when
                          considering samples with unequal coverage across
                          chromosomes, like male samples. An usage examples is
                          --ignoreForNormalization chrX chrM. (default: None)
    --skipNonCoveredRegions, --skipNAs
                          This parameter determines if non-covered regions
                          (regions without overlapping reads) in a BAM file
                          should be skipped. The default is to treat those
                          regions as having a value of zero. The decision to
                          skip non-covered regions depends on the interpretation
                          of the data. Non-covered regions may represent, for
                          example, repetitive regions that should be skipped.
                          (default: False)
    --smoothLength INT bp
                          The smooth length defines a window, larger than the
                          binSize, to average the number of reads. For example,
                          if the --binSize is set to 20 and the --smoothLength
                          is set to 60, then, for each bin, the average of the
                          bin and its left and right neighbors is considered.
                          Any value smaller than --binSize will be ignored and
                          no smoothing will be applied. (default: None)

  Read processing options:
    --extendReads [INT bp], -e [INT bp]
                          This parameter allows the extension of reads to
                          fragment size. If set, each read is extended, without
                          exception. *NOTE*: This feature is generally NOT
                          recommended for spliced-read data, such as RNA-seq, as
                          it would extend reads over skipped regions. *Single-
                          end*: Requires a user specified value for the final
                          fragment length. Reads that already exceed this
                          fragment length will not be extended. *Paired-end*:
                          Reads with mates are always extended to match the
                          fragment size defined by the two read mates. Unmated
                          reads, mate reads that map too far apart (>4x fragment
                          length) or even map to different chromosomes are
                          treated like single-end reads. The input of a fragment
                          length value is optional. If no value is specified, it
                          is estimated from the data (mean of the fragment size
                          of all mate reads). (default: False)
    --ignoreDuplicates    If set, reads that have the same orientation and start
                          position will be considered only once. If reads are
                          paired, the mate's position also has to coincide to
                          ignore a read. (default: False)
    --minMappingQuality INT
                          If set, only reads that have a mapping quality score
                          of at least this are considered. (default: None)
    --centerReads         By adding this option, reads are centered with respect
                          to the fragment length. For paired-end data, the read
                          is centered at the fragment length defined by the two
                          ends of the fragment. For single-end data, the given
                          fragment length is used. This option is useful to get
                          a sharper signal around enriched regions. (default:
                          False)
    --samFlagInclude INT  Include reads based on the SAM flag. For example, to
                          get only reads that are the first mate, use a flag of
                          64. This is useful to count properly paired reads only
                          once, as otherwise the second mate will be also
                          considered for the coverage. (default: None)
    --samFlagExclude INT  Exclude reads based on the SAM flag. For example, to
                          get only reads that map to the forward strand, use
                          --samFlagExclude 16, where 16 is the SAM flag for
                          reads that map to the reverse strand. (default: None)


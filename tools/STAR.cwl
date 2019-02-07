#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: InitialWorkDirRequirement
  listing:
    - writable: true
      entryname: generated_STAR_genome_dir
      entry: $(inputs.genomeDir)
hints:
- class: DockerRequirement
    #dockerImageId: scidap/star:v2.5.0b #not yet ready
  dockerPull: scidap/star:v2.5.0b
  dockerFile: >
    $import: STAR-Dockerfile

inputs:
  winBinNbits:
    type: int?
    inputBinding:
      position: 1
      prefix: --winBinNbits
    doc: '16

      int>0: =log2(winBin), where winBin is the size of the bin for the windows/clustering,
      each window will occupy an integer number of bins.

      '
  outFilterMatchNmin:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterMatchNmin
    doc: '0

      int: alignment will be output only if the number of matched bases is higher
      than this value

      '
  outSAMattributes:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMattributes
      shellQuote: false
    doc: |
      Standard
      string: a string of desired SAM attributes, in the order desired for the output SAM
      NH HI AS nM NM MD jM jI XS ... any combination in any order
      Standard   ... NH HI AS nM
      All        ... NH HI AS nM NM MD jM jI
      None       ... no attributes
  outSAMheaderPG:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --outSAMheaderPG
    doc: '-

      strings: extra @PG (software) line of the SAM header (in addition to STAR)

      '
  clip3pAdapterMMp:
    type: float?
    inputBinding:
      position: 1
      prefix: --clip3pAdapterMMp
    doc: |
      double(s): max proportion of mismatches for 3p adpater clipping for each
      mate.  If one value is given, it will be assumed the same for both mates.
  clip3pAfterAdapterNbases:
    type: int?
    inputBinding:
      position: 1
      prefix: --clip3pAfterAdapterNbases
    doc: |
      int(s): number of bases to clip from 3p of each mate after the adapter
      clipping. If one value is given, it will be assumed the same for both
      mates.
  scoreGapNoncan:
    type: int?
    inputBinding:
      position: 1
      prefix: --scoreGapNoncan
    doc: |
      -8
      int: non-canonical junction penalty (in addition to scoreGap)
  outMultimapperOrder:
    type: string?
    inputBinding:
      position: 1
      prefix: --outMultimapperOrder
    doc: |
      Old_2.4
      string: order of multimapping alignments in the output files
      Old_2.4             ... quasi-random order used before 2.5.0
      Random              ... random order of alignments for each multi-mapper. Read mates (pairs) are always adjacent, all alignment for each read stay together. This option will become default in the future releases.
  quantTranscriptomeBan:
    type: string?
    inputBinding:
      position: 1
      prefix: --quantTranscriptomeBan
    doc: |
      IndelSoftclipSingleend
      string: prohibit various alignment type
      IndelSoftclipSingleend  ... prohibit indels, soft clipping and single-end alignments - compatible with RSEM
      Singleend               ... prohibit single-end alignments
  alignSJstitchMismatchNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignSJstitchMismatchNmax
    doc: |
      0 -1 0 0
      4*int>=0: maximum number of mismatches for stitching of the splice junctions (-1: no limit).
      (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.
  twopassMode:
    type: string?
    inputBinding:
      position: 1
      prefix: --twopassMode
    doc: |
      None
      string: 2-pass mapping mode.
      None        ... 1-pass mapping
      Basic       ... basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly
  outSAMorder:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMorder
    doc: |
      Paired
      string: type of sorting for the SAM output
      Paired: one mate after the other for all paired alignments
      PairedKeepInputOrder: one mate after the other for all paired alignments, the order is kept the same as in the input FASTQ files
  outSJfilterDistToOtherSJmin:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --outSJfilterDistToOtherSJmin
    doc: |
      10  0   5   10
      4 integers>=0: minimum allowed distance to other junctions' donor/acceptor
      does not apply to annotated junctions
  parametersFiles:
    type: string?
    inputBinding:
      position: 1
      prefix: --parametersFiles
    doc: |
      string: name of a user-defined parameters file, "-": none. Can only be
      defined on the command line.
  limitSjdbInsertNsj:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitSjdbInsertNsj
    doc: '1000000

      int>=0: maximum number of junction to be inserted to the genome on the fly at
      the mapping stage, including those from annotations and those detected in the
      1st step of the 2-pass run

      '
  seedSearchLmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --seedSearchLmax
    doc: '0

      int>=0: defines the maximum length of the seeds, if =0 max seed lengthis infinite

      '
  bamRemoveDuplicatesMate2basesN:
    type: int?
    inputBinding:
      position: 1
      prefix: --bamRemoveDuplicatesMate2basesN
    doc: |
      int>0: number of bases from the 5' of mate 2 to use in collapsing (e.g. for
      RAMPAGE)
  alignTranscriptsPerWindowNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignTranscriptsPerWindowNmax
    doc: |
      100
      int>0: max number of transcripts per window
  outFilterIntronMotifs:
    type: string?
    inputBinding:
      position: 1
      prefix: --outFilterIntronMotifs
    doc: |
      None
      string: filter alignment using their motifs
      None                           ... no filtering
      RemoveNoncanonical             ... filter out alignments that contain non-canonical junctions
      RemoveNoncanonicalUnannotated  ... filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. The annotated non-canonical junctions will be kept.
  winAnchorDistNbins:
    type: int?
    inputBinding:
      position: 1
      prefix: --winAnchorDistNbins
    doc: '9

      int>0: max number of bins between two anchors that allows aggregation of anchors
      into one window

      '
  chimScoreDropMax:
    type: int?
    inputBinding:
      position: 1
      prefix: --chimScoreDropMax
    doc: '20

      int>=0: max drop (difference) of chimeric score (the sum of scores of all chimeric
      segements) from the read length

      '
  outSAMattrRGline:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMattrRGline
    doc: |
      -
      string(s): SAM/BAM read group line. The first word contains the read group identifier and must start with "ID:", e.g. --outSAMattrRGline ID:xxx CN:yy "DS:z z z".
      xxx will be added as RG tag to each output alignment. Any spaces in the tag values have to be double quoted.
      Comma separated RG lines correspons to different (comma separated) input files in --readFilesIn. Commas have to be surrounded by spaces, e.g.
      --outSAMattrRGline ID:xxx , ID:zzz "DS:z z" , ID:yyy DS:yyyy
  chimOutType:
    type: string?
    inputBinding:
      position: 1
      prefix: --chimOutType
    doc: |
      SeparateSAMold
      string: type of chimeric output
      SeparateSAMold  ... output old SAM into separate Chimeric.out.sam file
      WithinBAM       ... output into main aligned BAM files (Aligned.*.bam)
  runDirPerm:
    type: string?
    inputBinding:
      position: 1
      prefix: --runDirPerm
    doc: |
      User_RWX
      string: permissions for the directories created at the run-time.
      User_RWX ... user-read/write/execute
      All_RWX  ... all-read/write/execute (same as chmod 777)
  outQSconversionAdd:
    type: int?
    inputBinding:
      position: 1
      prefix: --outQSconversionAdd
    doc: |
      int: add this number to the quality score (e.g. to convert from Illumina to Sanger, use -31)
  outSAMattrIHstart:
    type: int?
    inputBinding:
      position: 1
      prefix: --outSAMattrIHstart
    doc: '1

      int>=0:                     start value for the IH attribute. 0 may be required
      by some downstream software, such as Cufflinks or StringTie.

      '
  chimSegmentMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --chimSegmentMin
    doc: '0

      int>=0: minimum length of chimeric segment length, if ==0, no chimeric output

      '
  scoreGapATAC:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --scoreGapATAC
    doc: |
      -8
      AT/AC  and GT/AT junction penalty  (in addition to scoreGap)
  readMatesLengthsIn:
    type: string?
    inputBinding:
      position: 1
      prefix: --readMatesLengthsIn
    doc: |
      string: Equal/NotEqual - lengths of names,sequences,qualities for both
      mates are the same  / not the same. NotEqual is safe in all situations.
  seedPerReadNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --seedPerReadNmax
    doc: |
      1000
      int>0: max number of seeds per read
  outWigType:
    type: string?
    inputBinding:
      position: 1
      prefix: --outWigType
    doc: |
      None
      string(s): type of signal output, e.g. "bedGraph" OR "bedGraph read1_5p". Requires sorted BAM: --outSAMtype BAM SortedByCoordinate .
      1st word:
      None       ... no signal output
      bedGraph   ... bedGraph format
      wiggle     ... wiggle format
      2nd word:
      read1_5p   ... signal from only 5' of the 1st read, useful for CAGE/RAMPAGE etc
      read2      ... signal from only 2nd read
  winFlankNbins:
    type: int?
    inputBinding:
      position: 1
      prefix: --winFlankNbins
    doc: '4

      int>0: log2(winFlank), where win Flank is the size of the left and right flanking
      regions for each window

      '
  sjdbGTFfeatureExon:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFfeatureExon
    doc: |
      exon
      string: feature type in GTF file to be used as exons for building
      transcripts
  genomeDir:
    type: Directory
    inputBinding:
      valueFrom: |
        ${
            if (inputs.runMode == "genomeGenerate")
              return "generated_STAR_genome_dir";
            return self;
        }
      position: 1
      prefix: --genomeDir
    doc: |
      string: path to the directory where genome files are stored (if
      runMode!=generateGenome) or will be generated (if runMode==generateGenome)
  chimFilter:
    type: string?
    inputBinding:
      position: 1
      prefix: --chimFilter
    doc: |
      banGenomicN
      string(s): different filters for chimeric alignments
      None ... no filtering
      banGenomicN ... Ns are not allowed in the genome sequence around the chimeric junction
  outSAMunmapped:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMunmapped
    doc: |
      string: output of unmapped reads in the SAM format
      None   ... no output
      Within ... output unmapped reads within the main SAM file (i.e. Aligned.out.sam)
  seedSearchStartLmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --seedSearchStartLmax
    doc: '50

      int>0: defines the search start point through the read - the read is split into
      pieces no longer than this value

      '
  runMode:
    type: string
    default: alignReads
    inputBinding:
      position: 1
      prefix: --runMode
    doc: |
      string: type of the run:
      alignReads             ... map reads
      genomeGenerate         ... generate genome files
      inputAlignmentsFromBAM ... input alignments from BAM. Presently only works with --outWigType and --bamRemoveDuplicates.
  genomeFastaFiles:
    type: File[]?
    inputBinding:
      position: 1
      itemSeparator: ' '
      prefix: --genomeFastaFiles
    doc: |
      string(s): path(s) to the fasta files with genomic sequences for genome
      generation, separated by spaces. Only used if runMode==genomeGenerate.
      These files should be plain text FASTA files, they *cannot* be zipped.
  limitIObufferSize:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitIObufferSize
    doc: |
      150000000
      int>0: max available buffers size (bytes) for input/output, per thread
  sjdbScore:
    type: int?
    inputBinding:
      position: 1
      prefix: --sjdbScore
    doc: |
      2
      int: extra alignment score for alignmets that cross database junctions
  alignSJDBoverhangMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignSJDBoverhangMin
    doc: '3

      int>0: minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments

      '
  outSAMstrandField:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMstrandField
    doc: |
      None
      string: Cufflinks-like strand field flag
      None        ... not used
      intronMotif ... strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.
  alignMatesGapMax:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignMatesGapMax
    doc: '0

      maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins

      '
  clip3pNbases:
    type: int?
    inputBinding:
      position: 1
      prefix: --clip3pNbases
    doc: |
      int(s): number(s) of bases to clip from 3p of each mate. If one value is
      given, it will be assumed the same for both mates.
  outFilterMultimapNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterMultimapNmax
    doc: '10

      int: read alignments will be output only if the read maps fewer than this value,
      otherwise no alignments will be output

      '
  outFileNamePrefix:
    type: string?
    inputBinding:
      position: 1
      prefix: --outFileNamePrefix
    doc: |
      string: output files name prefix (including full or relative path). Can
      only be defined on the command line.
  quantMode:
    type: string?
    inputBinding:
      position: 1
      prefix: --quantMode
    doc: |
      -
      string(s): types of quantification requested
      -                ... none
      TranscriptomeSAM ... output SAM/BAM alignments to transcriptome into a separate file
      GeneCounts       ... count reads per gene
  outFilterMismatchNoverLmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNoverLmax
    doc: '0.3

      int: alignment will be output only if its ratio of mismatches to *mapped* length
      is less than this value

      '
  sjdbGTFchrPrefix:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFchrPrefix
    doc: |
      string: prefix for chromosome names in a GTF file (e.g. 'chr' for using
      ENSMEBL annotations with UCSC geneomes)
  clip3pAdapterSeq:
    type: string?
    inputBinding:
      position: 1
      prefix: --clip3pAdapterSeq
    doc: |
      string(s): adapter sequences to clip from 3p of each mate.  If one value is
      given, it will be assumed the same for both mates.
  outBAMsortingThreadN:
    type: int?
    inputBinding:
      position: 1
      prefix: --outBAMsortingThreadN
    doc: |
      int: >=0: number of threads for BAM sorting. 0 will default to
      min(6,--runThreadN).
  twopass1readsN:
    type: int?
    inputBinding:
      position: 1
      prefix: --twopass1readsN
    doc: |
      int: number of reads to process for the 1st step. Use very large number (or
      default -1) to map all reads in the first step.
  outSJfilterIntronMaxVsReadN:
    type: int[]?
    inputBinding:
      position: 1
      prefix: --outSJfilterIntronMaxVsReadN
    doc: |
      50000 100000 200000
      N integers>=0: maximum gap allowed for junctions supported by 1,2,3,,,N reads
      i.e. by default junctions supported by 1 read can have gaps <=50000b, by 2 reads: <=100000b, by 3 reads: <=200000. by >=4 reads any gap <=alignIntronMax
      does not apply to annotated junctions
  outSAMfilter:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMfilter
    doc: |
      None
      string(s): filter the output into main SAM/BAM files
      KeepOnlyAddedReferences ... only keep the reads for which all alignments are to the extra reference sequences added with --genomeFastaFiles at the mapping stage.
  outSAMheaderHD:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --outSAMheaderHD
    doc: |
      -
      strings: @HD (header) line of the SAM header
  chimScoreMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --chimScoreMin
    doc: |
      0
      int>=0: minimum total (summed) score of the chimeric segments
  outSJfilterOverhangMin:
    type: int[]?
    inputBinding:
      position: 1
      prefix: --outSJfilterOverhangMin
    doc: |
      30  12  12  12
      4 integers:    minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif
      does not apply to annotated junctions
  scoreStitchSJshift:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --scoreStitchSJshift
    doc: '1

      maximum score reduction while searching for SJ boundaries inthe stitching step

      '
  readNameSeparator:
    type: string?
    inputBinding:
      position: 1
      prefix: --readNameSeparator
    doc: |
      /
      string(s): character(s) separating the part of the read names that will be
      trimmed in output (read name after space is always trimmed)
  scoreGapGCAG:
    type: int?
    inputBinding:
      position: 1
      prefix: --scoreGapGCAG
    doc: |
      -4
      GC/AG and CT/GC junction penalty (in addition to scoreGap)
  scoreInsBase:
    type: int?
    inputBinding:
      position: 1
      prefix: --scoreInsBase
    doc: |
      -2
      insertion extension penalty per base (in addition to scoreInsOpen)
  quantTranscriptomeBAMcompression:
    type: int?
    inputBinding:
      position: 1
      prefix: --quantTranscriptomeBAMcompression
    doc: |
      int: -1 to 10  transcriptome BAM compression level, -1=default compression
      (6?), 0=no compression, 10=maximum compression
  seedMultimapNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --seedMultimapNmax
    doc: '10000

      int>0: only pieces that map fewer than this value are utilized in the stitching
      procedure

      '
  genomeLoad:
    type: string?
    inputBinding:
      position: 1
      prefix: --genomeLoad
    doc: |
      NoSharedMemory
      string: mode of shared memory usage for the genome files
      LoadAndKeep     ... load genome into shared and keep it in memory after run
      LoadAndRemove   ... load genome into shared but remove it after run
      LoadAndExit     ... load genome into shared memory and exit, keeping the genome in memory for future runs
      Remove          ... do not map anything, just remove loaded genome from memory
      NoSharedMemory  ... do not use shared memory, each job will have its own private copy of the genome
  chimScoreJunctionNonGTAG:
    type: int?
    inputBinding:
      position: 1
      prefix: --chimScoreJunctionNonGTAG
    doc: |
      -1
      int: penalty for a non-GT/AG chimeric junction
  sjdbInsertSave:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbInsertSave
    doc: |
      Basic
      string: which files to save when sjdb junctions are inserted on the fly at the mapping step
      Basic ... only small junction / transcript files
      All   ... all files including big Genome, SA and SAindex - this will create a complete genome directory
  sjdbFileChrStartEnd:
    type: File[]?
    inputBinding:
      position: 1
      prefix: --sjdbFileChrStartEnd
    doc: '-

      string(s): path to the files with genomic coordinates (chr <tab> start <tab>
      end <tab> strand) for the splice junction introns. Multiple files can be supplied
      wand will be concatenated.

      '
  genomeChrBinNbits:
    type: int?
    inputBinding:
      position: 1
      prefix: --genomeChrBinNbits
    doc: |
      int: =log2(chrBin), where chrBin is the size of the bins for genome
      storage: each chromosome will occupy an integer number of bins
  readFilesIn:
    type: File[]?
    inputBinding:
      position: 1
      itemSeparator: ' '
      prefix: --readFilesIn
      shellQuote: false
    doc: |
      string(s): paths to files that contain input read1 (and, if needed,  read2)
  genomeSAsparseD:
    type: int?
    inputBinding:
      position: 1
      prefix: --genomeSAsparseD
    doc: |
      int>0: suffux array sparsity, i.e. distance between indices: use bigger
      numbers to decrease needed RAM at the cost of mapping speed reduction
  scoreDelOpen:
    type: int?
    inputBinding:
      position: 1
      prefix: --scoreDelOpen
    doc: |
      -2
      deletion open penalty
  outSAMflagAND:
    type: int?
    inputBinding:
      position: 1
      prefix: --outSAMflagAND
    doc: '65535

      int: 0 to 65535: sam FLAG will be bitwise AND''d with this value, i.e. FLAG=FLAG
      & outSAMflagOR. This is applied after all flags have been set by STAR, but before
      outSAMflagOR. Can be used to unset specific bits that are not set otherwise.

      '
  outSAMmultNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --outSAMmultNmax
    doc: '-1

      int: max number of multiple alignments for a read that will be output to the
      SAM/BAM files. -1 ... all alignments (up to --outFilterMultimapNmax) will be
      output

      '
  limitGenomeGenerateRAM:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitGenomeGenerateRAM
    doc: |
      31000000000
      int>0: maximum available RAM (bytes) for genome generation
  clip5pNbases:
    type: int?
    inputBinding:
      position: 1
      prefix: --clip5pNbases
    doc: |
      int(s): number(s) of bases to clip from 5p of each mate. If one value is
      given, it will be assumed the same for both mates.
  outWigReferencesPrefix:
    type: string?
    inputBinding:
      position: 1
      prefix: --outWigReferencesPrefix
    doc: |
      string: prefix matching reference names to include in the output wiggle
      file, e.g. "chr", default "-" - include all references
  alignTranscriptsPerReadNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignTranscriptsPerReadNmax
    doc: |
      10000
      int>0: max number of different alignments per read to consider
  outFilterScoreMinOverLread:
    type: float?
    inputBinding:
      position: 1
      prefix: --outFilterScoreMinOverLread
    doc: '0.66

      float: outFilterScoreMin normalized to read length (sum of mates'' lengths for
      paired-end reads)

      '
  outSAMreadID:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMreadID
    doc: |
      Standard
      string: read ID record type
      Standard ... first word (until space) from the FASTx read ID line, removing /1,/2 from the end
      Number   ... read number (index) in the FASTx file
  scoreDelBase:
    type: int?
    inputBinding:
      position: 1
      prefix: --scoreDelBase
    doc: |
      -2
      deletion extension penalty per base (in addition to scoreDelOpen)
  readMapNumber:
    type: int?
    inputBinding:
      position: 1
      prefix: --readMapNumber
    doc: |
      -1
      int: number of reads to map from the beginning of the file
      -1: map all reads
  inputBAMfile:
    type: File?
    inputBinding:
      position: 1
      prefix: --inputBAMfile
    doc: |
      string: path to BAM input file, to be used with --runMode
      inputAlignmentsFromBAM
  outSAMflagOR:
    type: int?
    inputBinding:
      position: 1
      prefix: --outSAMflagOR
    doc: |
      int: 0 to 65535: sam FLAG will be bitwise OR'd with this value, i.e.
      FLAG=FLAG | outSAMflagOR. This is applied after all flags have been set by
      STAR, and after outSAMflagAND. Can be used to set specific bits that are
      not set otherwise.
  sjdbGTFtagExonParentTranscript:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFtagExonParentTranscript
    doc: |
      transcript_id
      string: tag name to be used as exons' transcript-parents (default
      "transcript_id" works for GTF files)
  alignEndsType:
    type: string?
    inputBinding:
      position: 1
      prefix: --alignEndsType
    doc: |
      Local
      string: type of read ends alignment
      Local           ... standard local alignment with soft-clipping allowed
      EndToEnd        ... force end-to-end read alignment, do not soft-clip
      Extend5pOfRead1 ... fully extend only the 5p of the read1, all other ends: local alignment
  sjdbGTFtagExonParentGene:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFtagExonParentGene
    doc: 'gene_id

      string: tag name to be used as exons'' gene-parents (default "gene_id" works
      for GTF files)

      '
  outFilterMatchNminOverLread:
    type: float?
    inputBinding:
      position: 1
      prefix: --outFilterMatchNminOverLread
    doc: '0.66

      float: outFilterMatchNmin normalized to read length (sum of mates'' lengths
      for paired-end reads)

      '
  alignSplicedMateMapLminOverLmate:
    type: float?
    inputBinding:
      position: 1
      prefix: --alignSplicedMateMapLminOverLmate
    doc: |
      0.66
      float>0: alignSplicedMateMapLmin normalized to mate length
  outSJfilterReads:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSJfilterReads
    doc: |
      All
      string: which reads to consider for collapsed splice junctions output
      All: all reads, unique- and multi-mappers
      Unique: uniquely mapping reads only
  scoreGap:
    type: int?
    inputBinding:
      position: 1
      prefix: --scoreGap
    doc: |
      0
      int: splice junction penalty (independent on intron motif)
  outWigNorm:
    type: string?
    inputBinding:
      position: 1
      prefix: --outWigNorm
    doc: |
      RPM
      string: type of normalization for the signal
      RPM    ... reads per million of mapped reads
      None   ... no normalization, "raw" counts
  runThreadN:
    type: int?
    inputBinding:
      position: 1
      prefix: --runThreadN
    doc: |
      1
      int: number of threads to run STAR
  outSAMprimaryFlag:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMprimaryFlag
    doc: |
      OneBestScore
      string: which alignments are considered primary - all others will be marked with 0x100 bit in the FLAG
      OneBestScore ... only one alignment with the best score is primary
      AllBestScore ... all alignments with the best score are primary
  alignSJoverhangMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignSJoverhangMin
    doc: |
      5
      int>0: minimum overhang (i.e. block size) for spliced alignments
  chimJunctionOverhangMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --chimJunctionOverhangMin
    doc: |
      20
      int>=0: minimum overhang for a chimeric junction
  bamRemoveDuplicatesType:
    type: string?
    inputBinding:
      position: 1
      prefix: --bamRemoveDuplicatesType
    doc: |
      -
      string: mark duplicates in the BAM file, for now only works with sorted BAM feeded with inputBAMfile
      -               ... no duplicate removal/marking
      UniqueIdentical ... mark all multimappers, and duplicate unique mappers. The coordinates, FLAG, CIGAR must be identical
  seedNoneLociPerWindow:
    type: int?
    inputBinding:
      position: 1
      prefix: --seedNoneLociPerWindow
    doc: "10 \nint>0: max number of one seed loci per window\n"
  sjdbGTFfile:
    type: File?
    inputBinding:
      position: 1
      prefix: --sjdbGTFfile
    doc: |
      string: path to the GTF file with annotations
  alignWindowsPerReadNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignWindowsPerReadNmax
    doc: |
      10000
      int>0: max number of windows per read
  sysShell:
    type: string?
    inputBinding:
      position: 1
      prefix: --sysShell
    doc: |
      string: path to the shell binary, preferrably bash, e.g. /bin/bash.
      - ... the default shell is executed, typically /bin/sh. This was reported to fail on some Ubuntu systems - then you need to specify path to bash.
  outSJfilterCountUniqueMin:
    type: int[]?
    inputBinding:
      position: 1
      prefix: --outSJfilterCountUniqueMin
    doc: "3   1   1   1 \n4 integers: minimum uniquely mapping read count per junction\
      \ for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC\
      \ motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif\nJunctions\
      \ are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin\
      \ conditions are satisfied\ndoes not apply to annotated junctions\n"
  outFilterMultimapScoreRange:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterMultimapScoreRange
    doc: |
      1
      int: the score range below the maximum score for multimapping alignments
  alignSoftClipAtReferenceEnds:
    type: string?
    inputBinding:
      position: 1
      prefix: --alignSoftClipAtReferenceEnds
    doc: |
      Yes
      string: allow the soft-clipping of the alignments past the end of the chromosomes
      Yes ... allow
      No  ... prohibit, useful for compatibility with Cufflinks
  outSAMmode:
    type: string
    default: Full
    inputBinding:
      position: 1
      prefix: --outSAMmode
    doc: |
      string: mode of SAM output
      None ... no SAM output
      Full ... full SAM output
      NoQS ... full SAM but without quality scores
  outFilterType:
    type: string?
    inputBinding:
      position: 1
      prefix: --outFilterType
    doc: |
      Normal
      string: type of filtering
      Normal  ... standard filtering using only current alignment
      BySJout ... keep only those reads that contain junctions that passed filtering into SJ.out.tab
  alignIntronMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignIntronMin
    doc: '21

      minimum intron size: genomic gap is considered intron if its length>=alignIntronMin,
      otherwise it is considered Deletion

      '
  sjdbOverhang:
    type: int?
    inputBinding:
      position: 1
      prefix: --sjdbOverhang
    doc: '100

      int>0: length of the donor/acceptor sequence on each side of the junctions,
      ideally = (mate_length - 1)

      '
  alignIntronMax:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignIntronMax
    doc: '0

      maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins

      '
  limitOutSJcollapsed:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitOutSJcollapsed
    doc: |
      1000000
      int>0: max number of collapsed junctions
  outSAMtype:
    type:
      type: array
      items: string
    default: [BAM, SortedByCoordinate]
    inputBinding:
      position: 1
      prefix: --outSAMtype
    doc: |
      strings: type of SAM/BAM output
      1st word:
      BAM  ... output BAM without sorting
      SAM  ... output SAM without sorting
      None ... no SAM/BAM output
      2nd, 3rd:
      Unsorted           ... standard unsorted
      SortedByCoordinate ... sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM.
  scoreInsOpen:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --scoreInsOpen
    doc: |
      -2
      insertion open penalty
  outSAMmapqUnique:
    type: int?
    inputBinding:
      position: 1
      prefix: --outSAMmapqUnique
    doc: |
      255
      int: 0 to 255: the MAPQ value for unique mappers
  genomeSAindexNbases:
    type: int?
    inputBinding:
      position: 1
      prefix: --genomeSAindexNbases
    doc: |
      int: length (bases) of the SA pre-indexing string. Typically between 10 and
      15. Longer strings will use much more memory, but allow faster searches.
  winAnchorMultimapNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --winAnchorMultimapNmax
    doc: |
      50
      int>0: max number of loci anchors are allowed to map to
  limitOutSAMoneReadBytes:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitOutSAMoneReadBytes
    doc: '100000

      int>0: max size of the SAM record for one read. Recommended value: >(2*(LengthMate1+LengthMate2+100)*outFilterMultimapNmax

      '
  seedPerWindowNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --seedPerWindowNmax
    doc: |
      50
      int>0: max number of seeds per window
  readFilesCommand:
    type: string?
    inputBinding:
      position: 1
      prefix: --readFilesCommand
    doc: |
      string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout
      For example: zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc.
  outStd:
    type: string
    default: Log
    inputBinding:
      position: 1
      prefix: --outStd
    doc: |
      Log
      string: which output will be directed to stdout (standard out)
      Log                    ... log messages
      SAM                    ... alignments in SAM format (which normally are output to Aligned.out.sam file), normal standard output will go into Log.std.out
      BAM_Unsorted           ... alignments in BAM format, unsorted. Requires --outSAMtype BAM Unsorted
      BAM_SortedByCoordinate ... alignments in BAM format, unsorted. Requires --outSAMtype BAM SortedByCoordinate
      BAM_Quant              ... alignments to transcriptome in BAM format, unsorted. Requires --quantMode TranscriptomeSAM
  outSAMheaderCommentFile:
    type: string?
    inputBinding:
      position: 1
      prefix: --outSAMheaderCommentFile
    doc: |
      -
      string: path to the file with @CO (comment) lines of the SAM header
  outFilterMismatchNoverReadLmax:
    type: float?
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNoverReadLmax
    doc: '1

      int: alignment will be output only if its ratio of mismatches to *read* length
      is less than this value

      '
  outWigStrand:
    type: string?
    inputBinding:
      position: 1
      prefix: --outWigStrand
    doc: |
      Stranded
      string: strandedness of wiggle/bedGraph output
      Stranded   ...  separate strands, str1 and str2
      Unstranded ...  collapsed strands
  outSJfilterCountTotalMin:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --outSJfilterCountTotalMin
    doc: "3   1   1   1 \n4 integers: minimum total (multi-mapping+unique) read count\
      \ per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3)\
      \ GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that\
      \ motif\nJunctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin\
      \ conditions are satisfied\ndoes not apply to annotated junctions\n"
  scoreGenomicLengthLog2scale:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --scoreGenomicLengthLog2scale
    doc: '-0.25

      extra score logarithmically scaled with genomic length of the alignment: scoreGenomicLengthLog2scale*log2(genomicLength)

      '
  limitOutSJoneRead:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitOutSJoneRead
    doc: '1000

      int>0: max number of junctions for one read (including all multi-mappers)

      '
  runRNGseed:
    type: int?
    inputBinding:
      position: 1
      prefix: --runRNGseed
    doc: |
      777
      int: random number generator seed.
  outTmpDir:
    type: string?
    inputBinding:
      position: 1
      prefix: --outTmpDir
    doc: |
      string: path to a directory that will be used as temporary by STAR. All contents of this directory will be removed!
      - the temp directory will default to outFileNamePrefix_STARtmp
  outFilterScoreMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterScoreMin
    doc: '0

      int: alignment will be output only if its score is higher than this value

      '
  outReadsUnmapped:
    type: string?
    inputBinding:
      position: 1
      prefix: --outReadsUnmapped
    doc: |
      None
      string: output of unmapped reads (besides SAM)
      None    ... no output
      Fastx   ... output in separate fasta/fastq files, Unmapped.out.mate1/2
  outBAMcompression:
    type: int
    default: 10
    inputBinding:
      position: 1
      prefix: --outBAMcompression
    doc: |
      int: -1 to 10  BAM compression level, -1=default compression (6?), 0=no
      compression, 10=maximum compression
  limitBAMsortRAM:
    type: long?
    inputBinding:
      position: 1
      prefix: --limitBAMsortRAM
    doc: |
      int>=0: maximum available RAM for sorting BAM. If =0, it will be set to the
      genome index size. 0 value can only be used with --genomeLoad
      NoSharedMemory option.
  seedSearchStartLmaxOverLread:
    type: float?
    inputBinding:
      position: 1
      prefix: --seedSearchStartLmaxOverLread
    doc: '1.0

      float: seedSearchStartLmax normalized to read length (sum of mates'' lengths
      for paired-end reads)

      '
  chimScoreSeparation:
    type: int?
    inputBinding:
      position: 1
      prefix: --chimScoreSeparation
    doc: |
      int>=0: minimum difference (separation) between the best chimeric score and
      the next one
  outFilterMismatchNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNmax
    doc: '10

      int: alignment will be output only if it has fewer mismatches than this value

      '
  chimSegmentReadGapMax:
    type: int?
    inputBinding:
      position: 1
      prefix: --chimSegmentReadGapMax
    doc: |
      0
      int>=0: maximum gap in the read sequence between chimeric segments
  alignSplicedMateMapLmin:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignSplicedMateMapLmin
    doc: |
      0
      int>0: minimum mapped length for a read mate that is spliced
outputs:
  indices:
    type: Directory?
    outputBinding:
      glob: |
        ${
          if (inputs.runMode != "genomeGenerate")
            return [];
          return "generated_STAR_genome_dir";
        }

  aligned:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.runMode == "genomeGenerate")
            return [];

          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          if (inputs.outSAMtype.indexOf("SAM") > -1) {
              return p+"Aligned.out.sam";
          } else {
           if ( inputs.outSAMtype.indexOf("SortedByCoordinate") > -1 )
              return p+"Aligned.sortedByCoord.out.bam";
            else
              return p+"Aligned.out.bam";
          }
        }
    secondaryFiles: |
      ${
         var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
         return [
           {"path": p+"Log.final.out", "class":"File"},
           {"path": p+"SJ.out.tab", "class":"File"},
           {"path": p+"Log.out", "class":"File"}
         ];
      }

  mappingstats:
    type: File?
    outputBinding:
      loadContents: true
      glob: |
        ${
          if (inputs.runMode == "genomeGenerate")
            return [];

          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Log.final.out";
        }

  readspergene:
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"ReadsPerGene.out.tab";
        }

  transcriptomesam:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.quantMode != "TranscriptomeSAM")
            return null;
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Aligned.toTranscriptome.out.bam";
        }

  bamRemDups:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.bamRemoveDuplicatesType != "UniqueIdentical")
            return null;
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Processed.out.bam";
        }

baseCommand: [STAR]
$namespaces:
  s: http://schema.org/
$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  class: s:SoftwareSourceCode
  s:name: STAR
  s:about: 'Aligns RNA-seq reads to a reference genome using uncompressed suffix arrays.
    STAR has a potential for accurately aligning long (several kilobases) reads that
    are emerging from the third-generation sequencing technologies.

    '
  s:url: https://github.com/alexdobin/STAR
  s:codeRepository: https://github.com/alexdobin/STAR.git

  s:license:
  - https://opensource.org/licenses/GPL-3.0

  s:targetProduct:
    class: s:SoftwareApplication
    s:softwareVersion: 2.5.0b
    s:applicationCategory: commandline tool
  s:programmingLanguage: C++
  s:publication:
  - class: s:ScholarlyArticle
    id: https://doi.org/10.1093/bioinformatics/bts635

  s:author:
  - class: s:Person
    id: mailto:dobin@cshl.edu
    s:name: Alexander Dobin
    s:email: mailto:dobin@cshl.edu
#    foaf:fundedBy: "NHGRI (NIH) grant U54HG004557"
    s:worksFor:
    - class: s:Organization
      s:name: Cold Spring Harbor Laboratory, Cold Spring Harbor, NY, USA
s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/STAR.cwl
s:codeRepository: https://github.com/common-workflow-language/workflows

s:author:
  class: s:Person
  s:name: Andrey Kartashov
  s:email: mailto:Andrey.Kartashov@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab



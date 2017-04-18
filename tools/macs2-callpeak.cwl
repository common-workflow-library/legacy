#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: macs2-docker.yml
- class: InlineJavascriptRequirement

inputs:
# ---------------------------------- #
# ------ Input files arguments ----- #
# ---------------------------------- #
  control:
    type: File?
    inputBinding:
      position: 2
      prefix: --control
    doc: Control sample file.
  extsize:
    type: float?
    inputBinding:
      position: 1
      prefix: --extsize
    doc: 'The arbitrary extension size in bp. When nomodel is  true, MACS will use
      this value as fragment size to  extend each read towards 3'' end, then pile
      them up.  It''s exactly twice the number of obsolete SHIFTSIZE.  In previous
      language, each read is moved 5''->3''  direction to middle of fragment by 1/2
      d, then  extended to both direction with 1/2 d. This is  equivalent to say each
      read is extended towards 5''->3''  into a d size fragment. DEFAULT: 200. EXTSIZE
      and  SHIFT can be combined when necessary. Check SHIFT  option.'
  verbose:
    type: int?
    inputBinding:
      position: 1
      prefix: --verbose
    doc: |
      VERBOSE_LEVEL     Set verbose level of runtime message. 0: only show
      critical message, 1: show additional warning message,
      2: show process information, 3: show debug messages.
      DEFAULT:2
  seed:
    type: int?
    inputBinding:
      position: 1
      prefix: --seed
    doc: |
      SEED           Set the random seed while down sampling data. Must be
      a non-negative integer in order to be effective.
      DEFAULT: not set
  broad-cutoff:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --broad-cutoff
    doc: |
      BROADCUTOFF
      Cutoff for broad region. This option is not available
      unless --broad is set. If -p is set, this is a pvalue
      cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1
  bdg:
# ---------------------------------- #
# ------    Output arguments   ----- #
# ---------------------------------- #
    type: boolean?
    inputBinding:
      position: 1
      prefix: --bdg
    doc: '  Whether or not to save extended fragment pileup, and local lambda tracks
      (two files) at every bp into a bedGraph file. DEFAULT: True'
  ratio:
    type: float?
    inputBinding:
      position: 1
      prefix: --ratio
    doc: |
      RATIO         When set, use a custom scaling ratio of ChIP/control
      (e.g. calculated using NCIS) for linear scaling.
      DEFAULT: ingore
  call-summits:
# ------------------------------------- #
# ------ Shifting model arguments ----- #
# ------------------------------------- #
    type: boolean?
    inputBinding:
      position: 1
      prefix: --call-summits
    doc: 'If set, MACS will use a more sophisticated signal processing approach to
      find subpeak summits in each enriched peak region. DEFAULT: False '
  broad:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --broad
    doc: 'If set, MACS will try to call broad peaks by linking nearby highly enriched
      regions. The linking region is controlled by another cutoff through --linking-cutoff.
      The maximum linking region length is 4 times of d from MACS. DEFAULT: False '
  format:
    type: string?
    inputBinding:
      position: 1
      prefix: -f
    doc: '-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE}, --format
      {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE} Format of tag file,
      "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM"
      or "BOWTIE" or "BAMPE". The default AUTO option will let MACS decide which format
      the file is. Note that MACS can''t detect "BAMPE" or "BEDPE" format with "AUTO",
      and you have to implicitly specify the format for "BAMPE" and "BEDPE". DEFAULT:
      "AUTO".'
  buffer-size:
    type: int?
    inputBinding:
      position: 1
      prefix: --buffer-size
    doc: |
      BUFFER_SIZE
      Buffer size for incrementally increasing internal
      array size to store reads alignment information. In
      most cases, you don't have to change this parameter.
      However, if there are large number of
      chromosomes/contigs/scaffolds in your alignment, it's
      recommended to specify a smaller buffer size in order
      to decrease memory usage (but it will take longer time
      to read alignment files). Minimum memory requested for
      reading an alignment file is about # of CHROMOSOME *
      BUFFER_SIZE * 2 Bytes. DEFAULT: 100000
      Output arguments:
  nolambda:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --nolambda
    doc: 'If True, MACS will use fixed background lambda as local lambda for every
      peak region. Normally, MACS calculates a dynamic local lambda to reflect the
      local bias due to potential chromatin structure. '
  cutoff-analysis:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --cutoff-analysis
    doc: 'While set, MACS2 will analyze number or total length of peaks that can be
      called by different p-value cutoff then output a summary table to help user
      decide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt
      file. Note, minlen and maxgap may affect the results. WARNING: May take ~30
      folds longer time to finish. DEFAULT: False Post-processing options: '
  treatment:
    type:
      type: array
      items: File
    inputBinding:
      position: 2
      prefix: --treatment
    doc: 'Treatment sample file(s). If multiple files are given as -t A B C, then
      they will all be read and pooled together. IMPORTANT: the first sample will
      be used as the outputs basename.'
  slocal:
    type: int?
    inputBinding:
      position: 1
      prefix: --slocal
    doc: |
      SMALLLOCAL   The small nearby region in basepairs to calculate
      dynamic lambda. This is used to capture the bias near
      the peak summit region. Invalid if there is no control
      data. If you set this to 0, MACS will skip slocal
      lambda calculation. *Note* that MACS will always
      perform a d-size local lambda calculation. The final
      local bias should be the maximum of the lambda value
      from d, slocal, and llocal size windows. DEFAULT: 1000
  nomodel:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --nomodel
    doc: "\t Whether or not to build the shifting model. If True,  MACS will not build\
      \ model. by default it means  shifting size = 100, try to set extsize to change\
      \ it.  DEFAULT: False"
  keep-dup:
    type: string?
    inputBinding:
      position: 1
      prefix: --keep-dup
    doc: |
      KEEPDUPLICATES
      It controls the MACS behavior towards duplicate tags
      at the exact same location -- the same coordination
      and the same strand. The 'auto' option makes MACS
      calculate the maximum tags at the exact same location
      based on binomal distribution using 1e-5 as pvalue
      cutoff; and the 'all' option keeps every tags. If an
      integer is given, at most this number of tags will be
      kept at the same location. The default is to keep one
      tag at the same location. Default: 1
  down-sample:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --down-sample
    doc: 'When set, random sampling method will scale down the bigger sample. By default,
      MACS uses linear scaling. Warning: This option will make your result unstable
      and irreproducible since each time, random reads would be selected. Consider
      to use ''randsample'' script instead. <not implmented>If used together with
      --SPMR, 1 million unique reads will be randomly picked.</not implemented> Caution:
      due to the implementation, the final number of selected reads may not be as
      you expected! DEFAULT: False '
  trackline:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --trackline
    doc: 'Tells MACS to include trackline with bedGraph files. To include this trackline
      while displaying bedGraph at UCSC genome browser, can show name and description
      of the file as well. However my suggestion is to convert bedGraph to bigWig,
      then show the smaller and faster binary bigWig file at UCSC genome browser,
      as well as downstream analysis. Require --bdg to be set. Default: Not include
      trackline. '
  bw:
    type: int?
    inputBinding:
      position: 1
      prefix: --bw
    doc: |
      BW               Band width for picking regions to compute fragment
      size. This value is only used while building the
      shifting model. DEFAULT: 300
  fe-cutoff:
    type: float?
    inputBinding:
      position: 1
      prefix: --fe-cutoff
    doc: |
      FECUTOFF  When set, the value will be used to filter out peaks
      with low fold-enrichment. Note, MACS2 use 1.0 as
      pseudocount while calculating fold-enrichment.
      DEFAULT: 1.0
  fix-bimodal:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --fix-bimodal
    doc: 'Whether turn on the auto pair model process. If set, when MACS failed to
      build paired model, it will use the nomodel settings, the --exsize parameter
      to extend each tags towards 3'' direction. Not to use this automate fixation
      is a default behavior now. DEFAULT: False '
  to-large:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --to-large
    doc: 'When set, scale the small sample up to the bigger sample. By default, the
      bigger dataset will be scaled down towards the smaller dataset, which will lead
      to smaller p/qvalues and more specific results. Keep in mind that scaling down
      will bring down background noise more. DEFAULT: False '
  g:
    type: string?
    inputBinding:
      position: 1
      prefix: -g
    doc: Effective genome size. It can be 1.0e+9 or 1000000000,  or shortcuts:'hs'
      for human (2.7e9), 'mm' for mouse  (1.87e9), 'ce' for C. elegans (9e7) and 'dm'
      for  fruitfly (1.2e8), Default:hs.
  shift:
    type: int?
    inputBinding:
      position: 1
      prefix: --shift
    doc: '(NOT the legacy --shiftsize option!) The arbitrary shift in bp. Use discretion
      while setting it other than default value. When NOMODEL is set, MACS will use
      this value to move cutting ends (5'') towards 5''->3'' direction then apply
      EXTSIZE to extend them to fragments. When this value is negative, ends will
      be moved toward 3''->5'' direction. Recommended to keep it as default 0 for
      ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for
      detecting enriched cutting loci such as certain DNAseI-Seq datasets. Note, you
      can''t set values other than 0 if format is BAMPE for paired-end data. DEFAULT:
      0.'
  m:
    type: string?
    inputBinding:
      position: 1
      prefix: -m
    doc: |
      MFOLD MFOLD, --mfold MFOLD MFOLD
      Select the regions within MFOLD range of high-
      confidence enrichment ratio against background to
      build model. Fold-enrichment in regions must be lower
      than upper limit, and higher than the lower limit. Use
      as "-m 10 30". DEFAULT:5 50
  q:
# ------------------------------------- #
# ------ Shifting model arguments ----- #
# ------------------------------------- #
    type: float?
    inputBinding:
      position: 1
      prefix: -q
    doc: 'Minimum FDR (q-value) cutoff for peak detection. DEFAULT: 0.05. -q, and
      -p are mutually exclusive.'
  p:
    type: float?
    inputBinding:
      position: 1
      prefix: -p
    doc: 'Pvalue cutoff for peak detection. DEFAULT: not set.  -q, and -p are mutually
      exclusive. If pvalue cutoff is  set, qvalue will not be calculated and reported
      as -1  in the final .xls file..'
  s:
# ------------------------------------- #
# ------ Shifting model arguments ----- #
# ------------------------------------- #
    type: int?
    inputBinding:
      position: 1
      prefix: -s
    doc: |
      TSIZE, --tsize TSIZE
      Tag size. This will overide the auto detected tag
      size. DEFAULT: Not set
  llocal:
    type: int?
    inputBinding:
      position: 1
      prefix: --llocal
    doc: |
      LARGELOCAL   The large nearby region in basepairs to calculate
      dynamic lambda. This is used to capture the surround
      bias. If you set this to 0, MACS will skip llocal
      lambda calculation. *Note* that MACS will always
      perform a d-size local lambda calculation. The final
      local bias should be the maximum of the lambda value
      from d, slocal, and llocal size windows. DEFAULT:
      10000.
  SPMR:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --SPMR
    doc: 'If True, MACS will save signal per million reads for fragment pileup profiles.
      Require --bdg to be set. Default: False '
outputs:
  output_peak_file:
    type: File
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_peaks.*Peak')
      outputEval: $(self[0])
    doc: Peak calling output file in narrowPeak format.
  output_peak_xls_file:
    type: File
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_peaks.xls')
    doc: Peaks information/report file.
  output_peak_summits_file:
    type: File
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_summits.bed')

    doc: Peaks summits bedfile.
  output_ext_frag_bdg_file:
    type: File?
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_treat_pileup.bdg')
    doc: Bedgraph with extended fragment pileup.
baseCommand:
- macs2
- callpeak

arguments:
- valueFrom: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
    ''))
  prefix: -n
  position: 1


class: CommandLineTool
id: MaskPrimers.py
label: pRESTO MaskPrimers.py

description: |
  usage: MaskPrimers.py [-h] [-v]  ...
  
  Removes primers and annotates sequences with primer and barcode identifiers
  
  optional arguments:
    -h, --help     show this help message and exit
    -v, --version  show program's version number and exit
  
  subcommands:
                   Alignment method
      align        Find primer matches using pairwise local alignment
      score        Find primer matches by scoring primers at a fixed position
  
  output files:
    mask-pass      processed reads with successful primer matches.
    mask-fail      raw reads failing primer identification.
  
  output annotation fields:
    SEQORIENT      the orientation of the output sequence. Either F (input)
                   or RC (reverse complement of input).
    PRIMER         name of the best primer match.
    BARCODE        the sequence preceding the primer match. Only output when
                   the --barcode flag is specified.
  
requirements:

inputs:
  - id: "#subcommand"
    type: enum
    symbols: ["score", "align"]
    inputBinding:
      position: 0 

  - id: "#input_files"
    type: File 
    inputBinding:
      prefix: "-s"

  - id: "#primer_file"
    type: File 
    inputBinding:
      prefix: "-p"

  - id: "#mode"
    type: enum
    symbols: ["cut","mask","trim","tag"]
    inputBinding:
      prefix: "--mode"

  - id: "#barcode"
    type: boolean
    inputBinding:
      prefix: "--barcode"

  - id: "#start"
    type: int
    inputBinding:
      prefix: "--start"

  - id: "#maxerror"
    type: float
    inputBinding:
      prefix: "--maxerror"

  - id: "#nproc"
    type: int
    inputBinding:
      prefix: "--nproc"

  - id: "#outname"
    type: string
    inputBinding:
      prefix: "--outname"

  - id: "#outdir"
    type: string
    inputBinding:
      prefix: "--outdir"

  - id: "#pipeline_log"
    type: string
    inputBinding:
      prefix: "--outdir"

  - id: "#error_log"
    type: string
    inputBinding:
      prefix: "--outdir"


outputs:
  - id: "#outputs"
    type: 

# _primers-pass.fastq

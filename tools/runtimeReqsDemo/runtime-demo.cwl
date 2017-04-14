#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: ExpressionTool

requirements:
 - class: InlineJavascriptRequirement

inputs: []

outputs:
  outdir: string
  tmpdir: string
  cores: int
  ram: int
  outdirSize: int
  tmpdirSize: int

expression: >
  ${ return {
   "outdir": runtime.outdir,
   "tmpdir": runtime.tmpdir,
   "cores": runtime.cores,
   "ram": runtime.ram,
   "outdirSize": runtime.outdirSize,
   "tmpdirSize": runtime.tmpdirSize } }  

#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
 - class: ScatterFeatureRequirement
 - class: InlineJavascriptRequirement
 - class: SchemaDefRequirement
   types:
    - name: FilePairs
      type: record
      fields:
        - name: forward
          type: File
        - name: reverse
          type: File

inputs:
  fastqSeqs: FilePairs[]

outputs:
  reports:
    type: Directory[]
    outputSource: runFastqc/report
  files:
    type: File[]
    outputSource: fastqSeqs

steps:
  arrayOfFilePairsToFileArray:
    run:
      class: ExpressionTool
      inputs:
        arrayOfFilePairs: FilePairs[]
      outputs:
        pairByPairs: File[]
      expression: >
        ${
        var val;
        var ret = [];
        for (val of inputs.arrayOfFilePairs) {
          ret.push(val.forward);
          ret.push(val.reverse);
        }
        return { 'pairByPairs': ret } ; }
    in:
      arrayOfFilePairs: fastqSeqs
    out: [ pairByPairs ]

  runFastqc:
    run: ../../tools/fastqc.cwl
    in:
      fastqFile: arrayOfFilePairsToFileArray/pairByPairs
    scatter: fastqFile
    out: [ report ]

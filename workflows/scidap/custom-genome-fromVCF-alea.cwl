#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Dr. Barski Lab / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: Workflow

cwlVersion: v1.0

inputs:
- id: reference
  type: File
  doc: |
    the reference genome fasta file
- id: phasedsnps
  type: File
  doc: |
    the phased SNPs variants vcf file
- id: phasedindels
  type: File
  doc: |
    the phased Indels variants vcf file
- id: strain
  type: string

- id: filename
  type: string

outputs:

- id: outfile
  type: File
  outputSource: '#applyindels/output'

steps:

- id: applysnps
  run: ../../tools/alea-insilico.cwl
  inputs:
  - {id: reference, source: '#reference'}
  - {id: phased, source: '#phasedsnps'}
  - {id: strain, source: '#strain'}
  - {id: output_filename, default: intermediate.fasta}
  outputs:
  - {id: output}

- id: applyindels
  run: ../../tools/alea-insilico.cwl
  inputs:
  - {id: reference, source: '#applysnps/output'}
  - {id: phased, source: '#phasedindels'}
  - {id: strain, source: '#strain'}
  - {id: output_filename, source: '#filename'}
  outputs:
  - {id: output}
doc: creates custom genome from reference genome and two phased VCF files SNPs and
  Indels


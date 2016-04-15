#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: samtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#INTSTR[]'
    type: int
    description: 'INT|STR[,...]'
    inputBinding:
      position: 2
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
  - samtools
  - flags
description: |
  About: Convert between textual and numeric flag representation
  Usage: samtools flags INT|STR[,...]

  Flags:
  	0x1	PAIRED        .. paired-end (or multiple-segment) sequencing technology
  	0x2	PROPER_PAIR   .. each segment properly aligned according to the aligner
  	0x4	UNMAP         .. segment unmapped
  	0x8	MUNMAP        .. next segment in the template unmapped
  	0x10	REVERSE       .. SEQ is reverse complemented
  	0x20	MREVERSE      .. SEQ of the next segment in the template is reversed
  	0x40	READ1         .. the first segment in the template
  	0x80	READ2         .. the last segment in the template
  	0x100	SECONDARY     .. secondary alignment
  	0x200	QCFAIL        .. not passing quality controls
  	0x400	DUP           .. PCR or optical duplicate
  	0x800	SUPPLEMENTARY .. supplementary alignment


#!/usr/bin/env cwl-runner
#
# Auto generated CWL please use with caution
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center
# Developed for CWL consortium http://commonwl.org/

class: CommandLineTool
requirements:
  - import: node-engine.cwl
  - import: envvar-global.yml
  - import: bedtools-docker.yml
inputs:
  - id: '#stdoutfile'
    type: string
  - id: '#outhtml'
    type: File
    description: out.html
    inputBinding:
      position: 4
  - id: '#'
    type: boolean
    description: '>'
    inputBinding:
      position: 3
  - id: '#i'
    type: File
    description: '<bed/gff/vcf>'
    inputBinding:
      position: 2
      prefix: '-i'
  - id: '#base'
    type:
      - 'null'
      - boolean
    description: "The browser basename.  Default: http://genome.ucsc.edu \n"
    inputBinding:
      position: 1
      prefix: '-base'
  - id: '#org'
    type:
      - 'null'
      - boolean
    description: |
      The organism. Default: human
    inputBinding:
      position: 1
      prefix: '-org'
  - id: '#db'
    type:
      - 'null'
      - boolean
    description: |
      The build.  Default: hg18
      Example:
      By default, the links created will point to human (hg18) UCSC browser.
      If you have a local mirror, you can override this behavior by supplying
      the -base, -org, and -db options.
      For example, if the URL of your local mirror for mouse MM9 is called:
      http://mymirror.myuniversity.edu, then you would use the following:
    inputBinding:
      position: 1
      prefix: '-db'
  - id: '#base'
    type:
      - 'null'
      - boolean
    description: |
      http://mymirror.myuniversity.edu
    inputBinding:
      position: 1
      prefix: '-base'
  - id: '#org'
    type:
      - 'null'
      - boolean
    description: |
      mouse
    inputBinding:
      position: 1
      prefix: '-org'
  - id: '#db'
    type:
      - 'null'
      - boolean
    description: |
      mm9
    inputBinding:
      position: 1
      prefix: '-db'
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
  - bedtools
  - links
description: |
  Tool:    bedtools links (aka linksBed)
  Version: v2.25.0
  Summary: Creates HTML links to an UCSC Genome Browser from a feature file.

  Usage:   bedtools links [OPTIONS] -i <bed/gff/vcf> > out.html

  Options: 
  	-base	The browser basename.  Default: http://genome.ucsc.edu 
  	-org	The organism. Default: human
  	-db	The build.  Default: hg18

  Example: 
  	By default, the links created will point to human (hg18) UCSC browser.
  	If you have a local mirror, you can override this behavior by supplying
  	the -base, -org, and -db options.

  	For example, if the URL of your local mirror for mouse MM9 is called: 
  	http://mymirror.myuniversity.edu, then you would use the following:
  	-base http://mymirror.myuniversity.edu
  	-org mouse
  	-db mm9


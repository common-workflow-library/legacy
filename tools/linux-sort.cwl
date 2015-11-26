#!/usr/bin/env cwl-runner

adms:Asset:
  admssw:SoftwareProject:
    doap:name: "sort"
    doap:description: >
      sort - sort lines of text files
    doap:release:
      - doap:revision: "5.93"
    doap:license: "GPL"
    doap:category: "commandline tool"
    doap:developer:
      - foaf:Person:
        foaf:name: "Mike Haertel"
      - foaf:Person:
        foaf:name: "Paul Eggert"
      - foaf:Person:
        foaf:mbox: "bug-coreutils@gnu.org"
  adms:AssetDistribution:
    doap:name: "linux-sort.cwl"
    doap:description: "Developed for CWL consortium http://commonwl.org/"
    doap:specification: "http://common-workflow-language.github.io/draft-3/"
    doap:release: "cwl:draft-3.dev2"
    doap:homepage: "http://commonwl.org/"
    doap:location: "https://github.com/common-workflow-language/workflows/blob/master/tools/linux-sort.cwl"
    doap:repository:
      - doap:GitRepository:
        doap:location: "https://github.com/common-workflow-language/workflows"
    doap:maintainer:
      foaf:Person:
        foaf:openid: "http://orcid.org/0000-0001-9102-5681"
        foaf:name: "Andrey Kartashov"
        foaf:mbox: "mailto:Andrey.Kartashov@cchmc.org"
        foaf:organization: "Cincinnati Children's Hospital Medical Center"

cwlVersion: "cwl:draft-3.dev2"

class: CommandLineTool

description: |
  Usage: sort [OPTION]... [FILE]...
    or:  sort [OPTION]... --files0-from=F
  Write sorted concatenation of all FILE(s) to standard output.

requirements:
  - "@import": envvar-global.cwl
  - "@import": linux-sort-docker.cwl
  - class: InlineJavascriptRequirement

inputs:
  - id: "#input"
    type: File
    inputBinding:
      position: 4

  - id: "#key"
    type: 
      type: array
      items: string
      inputBinding:
        prefix: "-k"
    inputBinding:
      position: 1

    description: |
      -k, --key=POS1[,POS2]
      start a key at POS1, end it at POS2 (origin 1)

outputs:
  - id: "#sorted"
    type: File
    description: "The sorted file"
    outputBinding:
      glob: $(inputs.input.path.split('/').slice(-1)[0] + '.sorted')

stdout: $(inputs.input.path.split('/').slice(-1)[0] + '.sorted')

baseCommand: ["sort"]

#arguments:
#  - valueFrom:
#      engine: node-engine.cwl
#      script: $job['key'].map(function(i) {return "-k"+i;})


#!/usr/bin/env cwl-runner

class: CommandLineTool
id: BAMStats
label: BAMStats tool
cwlVersion: v1.0

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  adms: http://www.w3.org/ns/adms#
  dcat: http://www.w3.org/ns/dcat#

$schemas:
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf

doc: |
  ![build_status](https://quay.io/repository/briandoconnor/dockstore-tool-bamstats/status)
  A Docker container for the BAMStats command. See the [BAMStats](http://bamstats.sourceforge.net/) website for more information.
  ```
  Usage:
  # fetch CWL
  $> dockstore tool cwl --entry quay.io/briandoconnor/dockstore-tool-bamstats:1.25-7 > Dockstore.cwl
  # make a runtime JSON template and edit it (or use the content of sample_configs.json in this git repo)
  $> dockstore tool convert cwl2json --cwl Dockstore.cwl > Dockstore.json
  # run it locally with the Dockstore CLI
  $> dockstore tool launch --entry quay.io/briandoconnor/dockstore-tool-bamstats:1.25-7 \
      --json Dockstore.json
  ```

dct:creator:
  '@id': http://orcid.org/0000-0002-7681-6415
  foaf:name: Brian O'Connor
  foaf:mbox: mailto:briandoconnor@gmail.com

requirements:
- class: DockerRequirement
  dockerPull: quay.io/briandoconnor/dockstore-tool-bamstats:1.25-7
- class: InlineJavascriptRequirement

hints:
- class: ResourceRequirement
  coresMin: 1
  ramMin: 4092
  outdirMin: 512000
  description: the process requires at least 4G of RAM

inputs:
  bam_input:
    type: File
    format: http://edamontology.org/format_2572
    inputBinding:
      position: 2
    doc: The BAM file used as input, it must be sorted.

  mem_gb:
    type: int
    default: 4
    inputBinding:
      position: 1
    doc: The memory, in GB, for the reporting tool

outputs:
  bamstats_report:
    type: File
    format: http://edamontology.org/format_3615
    outputBinding:
      glob: bamstats_report.zip
    doc: A zip file that contains the HTML report and various graphics.

baseCommand: [bash, /usr/local/bin/bamstats]


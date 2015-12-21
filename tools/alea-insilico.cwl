#!/usr/bin/env cwl-runner
#
# Author: Andrey.Kartashov@cchmc.org (http://orcid.org/0000-0001-9102-5681) / Cincinnati Children’s Hospital Medical Center / Dr. Barski Lab
# Developed for CWL consortium http://commonwl.org/

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

description: |
  alea-insilico.cwl provides direct call of alea.jar with insilico param

requirements:
- $import: envvar-global.cwl
- $import: alea-docker.cwl
- class: InlineJavascriptRequirement

inputs:
- id: "reference"
  type: File
  description: |
    the reference genome fasta file
  inputBinding:
    separate: false
    prefix: "--input-fasta="
    position: 2
    secondaryFiles:
    - ".fai"

- id: "phased"
  type: File
  description: |
    the phased variants vcf file (SNPs or/and Indels)
  inputBinding:
    separate: false
    prefix: "--input-vcf="
    position: 3
    secondaryFiles:
    - ".tbi"

- id: "strain"
  type: string
  description: |
    name of strain1 exactly as specified in the vcf file (e.g. hap1)
    mouse:
        129P2_OlaHsd	(129P2/OlaHsd)	F	52
        129S1_SvImJ	(129S1/SvImJ)	F	68
        129S5SvEvBrd	(129S5SvEvBrd)	F	22
        A_J	(A/J)	F	52
        AKR_J	(AKR/J)	F	57
        BALB_cJ	(BALB/cJ)	F	62
        BTBR_T+_Itpr3tf_J	(BTBR T+ Itpr3tf/J)	M	85
        BUB_BnJ	(BUB/BnJ)	M	49
        C3H_HeH	(C3H/HeH)	F	14
        C3H_HeJ	(C3H/HeJ)	F	63
        C57BL_10J	(C57BL/10J)	M	37
        C57BL_6NJ	(C57BL/6NJ)	F	61
        C57BR_cdJ	(C57BR/cdJ)	M	51
        C57L_J	(C57L/J)	M	64
        C58_J	(C58/J)	M	55
        CAST_EiJ	(CAST/EiJ)	F	53
        CBA_J	(CBA/J)	F	56
        DBA_1J	(DBA/1J)	M	49
        DBA_2J	(DBA/2J)	F	56
        FVB_NJ	(FVB/NJ)	F	73
        I_LnJ	(I/LnJ)	M	45
        KK_HiJ	(KK/HiJ)	M	55
        LEWES_EiJ	(LEWES/EiJ)	F	19
        LP_J	(LP/J)	F	54
        MOLF_EiJ	(MOLF/EiJ)	M	40
        NOD_ShiLtJ	(NOD/ShiLtJ)	F	66
        NZB_B1NJ	(NZB/B1NJ)	M	47
        NZO_HlLtJ	(NZO/HlLtJ)	F	72
        NZW_LacJ	(NZW/LacJ)	M	58
        PWK_PhJ	(PWK/PhJ)	F	53
        RF_J	(RF/J)	M	54
        SEA_GnJ	(SEA/GnJ)	M	49
        SPRET_EiJ	(SPRET/EiJ)	F	67
        ST_bJ	(ST/bJ)	M	81
        WSB_EiJ	(WSB/EiJ)	F	51
        ZALENDE_EiJ	(ZALENDE/EiJ)	M	19
  inputBinding:
    separate: false
    prefix: "--strain="
    position: 4

- id: "output_filename"
  type: string
  inputBinding:
    position: 5
    separate: false
    prefix: "--output-fasta="

outputs:
- id: "output"
  type: File
  outputBinding:
    glob: $(inputs.output_filename)
    secondaryFiles:
    - ".fai"
    - ".refmap"

baseCommand: ["java", "-Xms4G", "-Xmx8G", "-jar", "/usr/local/bin/alea.jar" ,"insilico"]

$namespaces:
  schema: http://schema.org/
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  adms: http://www.w3.org/ns/adms#
  dcat: http://www.w3.org/ns/dcat#

$schemas:
- http://schema.rdfs.org/all.rdf
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf

adms:includedAsset:
  $include: alea-ontology.yaml

dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/alea-insilico.cwl"
doap:repository:
- class: doap:GitRepository
  doap:location: "https://github.com/common-workflow-language/workflows"
doap:homepage: "http://commonwl.org/"
doap:license: "Apache2"


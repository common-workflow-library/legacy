#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool


doc: |
      GATK-VariantsRecalibrator.cwl is developed for CWL consortium

      Usage:
      ```
      java -Xmx8G \
            -jar gatk.jar
            -T VariantRecalibrator \
            -R [reference_fasta] \
            -recalFile    $tmpDir/out.recal \
            -tranchesFile $tmpDir/out.tranches \
            -rscriptFile  $tmpDir/out.R \
            -nt 4 \
            -an MQRankSum -an ReadPosRankSum -an DP -an FS -an QD \
            -mode SNP \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 [hapmap_vcf] \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0  [dbsnp_vcf] \
            -resource:omni,known=false,training=true,truth=true,prior=12.0   [1komni_vcf] \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 [1ksnp_vcf]
      ```


hints:
- $import: envvar-global.yml
- $import: GATK-docker.yml

inputs:
  haplotypecaller_snps_vcf:
    type: File
    inputBinding:
      position: 5
      prefix: -input
    doc: input vcf File raw variant calls from haplotype caller

  multithreading_nt:
    type: int
    default: 1
    inputBinding:
      position: 6
      prefix: -nt
    doc: multithreading option

  reference:
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    inputBinding:
      position: 7
      prefix: -R
    doc: reference genome

  resource_hapmap:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 8
      prefix: "-resource:hapmap,known=false,training=true,truth=true,prior=15.0"
    doc: hapmap reference data

  resource_omni:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 9
      prefix: "-resource:omni,known=false,training=true,truth=false,prior=12.0"
    doc: omni reference data

  resource_mills:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 9
      prefix: "-resource:mills,known=false,training=true,truth=true,prior=12.0"
    doc: hapmap reference data

  resource_dbsnp:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 9
      prefix: "-resource:dbsnp,known=true,training=false,truth=false,prior=8.0"
    doc: dbSNP reference data

  max_gaussian:
    type: int
    default: 4
    inputBinding:
      position: 10
      prefix: --maxGaussians

  java_arg:
    type: string
    default: -Xmx8g
    inputBinding:
      position: 1

  java_tmp:
    type: string
    default: -Djava.io.tmpdir=/tmp
    inputBinding:
      position: 2

outputs:
    tranches_File:
      type: File
      outputBinding:
        glob: vqsr_tranches.indels.recal
      doc: the tranches File

    recal_File:
      type: File
      outputBinding:
        glob: vqsr_recal.indels.recal
      doc: the recal File

   # vqsr_rscript:
   #     type: File
   #     outputBinding:
   #       glob: vqsr_tranches.plots.R
   #     doc: The output recalibration R script for the plots


arguments:
#- valueFrom: ./test/test-Files
#  position: 2
#  separate: false
#  prefix: -Djava.io.tmpdir=
- valueFrom: /usr/local/bin/GenomeAnalysisTK.jar
  position: 3
  prefix: -jar

- valueFrom: VariantRecalibrator
  position: 4
  prefix: -T

- valueFrom: "INDEL"
  position: 11
  prefix: -mode

- valueFrom: "QD"
  position: 12
  prefix: -an

#- valueFrom: "MQ"
#  position: 12
#  prefix: -an

- valueFrom: "MQRankSum"
  position: 12
  prefix: -an

- valueFrom: "ReadPosRankSum"
  position: 12
  prefix: -an

- valueFrom: "FS"
  position: 12
  prefix: -an
#- valueFrom: "SOR"
#  position: 12
#  prefix: -an

- valueFrom: vqsr_tranches.indels.recal
  position: 13
  prefix: -tranchesFile
- valueFrom: vqsr_recal.indels.recal
  position: 14
  prefix: -recalFile
#- valueFrom: vqsr_tranches.indels.plots.R
#  position: 15
#  prefix: -rscriptFile

baseCommand: [java]

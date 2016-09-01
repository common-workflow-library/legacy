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

inputs: # position 0, for java args, 1 for the jar, 2 for the tool itself
  GATKJar:
    type: File
    inputBinding:
      position: 1
      prefix: "-jar"
  haplotypecaller_snps_vcf:
    type: File
#    secondaryFiles:
#      - .idx
    inputBinding:
      position: 2
      prefix: -input
    doc: input vcf File raw variant calls from haplotype caller

  multithreading_nt:
    type: int
    default: 4
    inputBinding:
      position: 2
      prefix: --nt
    doc: multithreading option

  reference:
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    inputBinding:
      position: 2
      prefix: -R
    doc: reference genome

  resource_hapmap:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 2
      prefix: "-resource:hapmap,known=false,training=true,truth=true,prior=15.0"
    doc: hapmap reference data

  resource_omni:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 2
      prefix: "-resource:omni,known=false,training=true,truth=false,prior=12.0"
    doc: omni reference data

  resource_1kg:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 2
      prefix: "-resource:1000G,known=false,training=true,truth=false,prior=10.0"
    doc: 1000 genome reference data

  resource_dbsnp:
    type: File
    secondaryFiles:
      - .idx
    inputBinding:
      position: 2
      prefix: "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0"
    doc: dbSNP  reference data


  java_arg:
    type: string
    default: -Xmx8g
    inputBinding:
      position: 0


outputs:
    tranches_File:
      type: File
      outputBinding:
        glob: vqsr_tranches.out
      doc: the tranches File

    recal_File:
      type: File
      outputBinding:
        glob: vqsr_recal.out
      doc: the recal File

    vqsr_rscript:
        type: File
        outputBinding:
          glob: vqsr.R
        doc: The output recalibration R script for the plots


arguments:
- valueFrom: $(runtime.tmpdir)
  position: 0
  separate: false
  prefix: -Djava.io.tmpdir=
- valueFrom: VariantRecalibrator
  position: 2
  prefix: -T
- valueFrom: "SNP"
  position: 2
  prefix: -mode
- valueFrom: "QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff"
  position: 2
  prefix: -an
- valueFrom: vqsr_tranches.out
  position: 2
  prefix: -tranchesFile
- valueFrom: vqsr_recal.out
  position: 2
  prefix: -recalFile
- valueFrom: vqsr.R
  position: 2
  prefix: -rscriptFile

baseCommand: [java]

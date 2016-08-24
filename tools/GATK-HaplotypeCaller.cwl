#!/usr/bin/env cwl-runner

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

cwlVersion: "cwl:draft-3"

class: CommandLineTool

adms:includedAsset:
  doap:name: "GATK"
  doap:description: >
    The Genome Analysis Toolkit or GATK is a software package for analysis of high-throughput sequencing data, developed by the Data Science and Data Engineering group at the Broad Institute. 
    The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance.
    Its robust architecture, powerful processing engine and high-performance computing features make it capable of taking on projects of any size.
    http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles
  doap:homepage: "https://www.broadinstitute.org/gatk/"
  doap:repository:
  - class: doap:GitRepository
    doap:location: "https://github.com/broadgsa/gatk.git"
  doap:release:
  - class: doap:Version
    doap:revision: "3.4"
  doap:license: "mixed licensing model"
  doap:category: "commandline tool"
  doap:programming-language: "JAVA"
  doap:developer:
  - class: foaf:Organization
    foaf:name: "Broad Institute"

description: |
  GATK-RealignTargetCreator.cwl is developed for CWL consortium
  Call germline SNPs and indels via local re-assembly of haplotypes

doap:name: "GATK-IndelRealigner.cwl"
dcat:downloadURL: "https://github.com/common-workflow-language/workflows/blob/master/tools/GATK-IndelRealigner.cwl"

dct:creator:
- class: foaf:Organization
  foaf:name: "THE UNIVERSITY OF MELBOURNE"
  foaf:member:
  - class: foaf:Person
    id: "farahk@student.unimelb.edu.au"
    foaf:name: "Farah Zaib Khan"
    foaf:mbox: "mailto:farahk@student.unimelb.edu.au"
    
  - class: foaf:Person
    id: "skanwal@student.unimelb.edu.au"
    foaf:name: "Sehrish Kanwal"
    foaf:mbox: "mailto:skanwal@student.unimelb.edu.au"

doap:maintainer:
- class: foaf:Organization
  foaf:name: "THE UNIVERSITY OF MELBOURNE"
  foaf:member:
  - class: foaf:Person
    id: "farahk@student.unimelb.edu.au"
    foaf:name: "Farah Zaib Khan"
    foaf:mbox: "mailto:farahk@student.unimelb.edu.au"
  - class: foaf:Person
    id: "skanwal@student.unimelb.edu.au"
    foaf:name: "Sehrish Kanwal"
    foaf:mbox: "mailto:skanwal@student.unimelb.edu.au"
    
requirements:
- $import: envvar-global.yml
- $import: envvar-global.yml
- $import: GATK-docker.yml

    
inputs:

  - id: "#java_arg"
    type: string
    default: "-Xmx4g"
    inputBinding: 
      position: 1
     
  - id: "#reference"
    type: File
    inputBinding:
      position: 5
      prefix: "-R"
    secondaryFiles:
      - ".amb"
      - ".ann"
      - ".bwt"
      - ".pac"
      - ".rbwt"
      - ".rpac"
      - ".rsa"
      - ".sa"
      - ".fai"
      - "^.dict"
        
  - id: "#inputBam_HaplotypeCaller"
    type: File
    description: bam file produced after printReads
    inputBinding:
      position: 6 
      prefix: "-I"
    secondaryFiles:
      - "^.bai"
    
  - id: "#dbsnp"  
    type: File
    description: latest_dbsnp.vcf set of known indels
    inputBinding: 
      position: 7
      prefix: "--dbsnp"
      
  - id: "#outputfile_HaplotypeCaller"
    type: ["null", string]
    description: name of the output file from HaplotypeCaller
    inputBinding:
      position: 8
      prefix: "-o"

  - id: "#useFilteredReadsForAnnotations"
    type: ["null", boolean]
    description: Use the contamination-filtered read maps for the purposes of annotating variants
    inputBinding:
      position: 9
      prefix: "--useFilteredReadsForAnnotations"

  - id: "#useAllelesTrigger"
    type: ["null", boolean]
    description: Use additional trigger on variants found in an external alleles file
    inputBinding:
      position: 10
      prefix: "--useAllelesTrigger"

  - id: "#stand_emit_conf"
    type: ["null", double]
    description: The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)
    inputBinding:
      position: 11
      prefix: "--standard_min_confidence_threshold_for_emitting"


  - id: "#stand_call_conf"
    type: ["null", double]
    description: The minimum phred-scaled confidence threshold at which variants should be called
    inputBinding:
      position: 12
      prefix: "--standard_min_confidence_threshold_for_calling"

  - id: "#sample_ploidy"
    type: ["null", int]
    description: Use additional trigger on variants found in an external alleles file
    inputBinding:
      position: 13
      prefix: "--sample_ploidy"

  - id: "#sample_name"
    type: ["null", string]
    description: Use additional trigger on variants found in an external alleles file
    inputBinding:
      position: 14
      prefix: "--sample_name"

  - id: "#globalMAPQ"
    type: ["null", int]
    description: The global assumed mismapping rate for reads
    inputBinding:
      position: 15
      prefix: "--phredScaledGlobalReadMismappingRate"

  - id: "#pcr_indel_model"
    type: ["null", string]
    description: The PCR indel model to use
    inputBinding:
      position: 16
      prefix: "--pcr_indel_model"

  - id: "#output_mode"
    type: ["null", string]
    description: The PCR indel model to use
    inputBinding:
      position: 17
      prefix: "--output_mode"
    
  - id: "#numPruningSamples"
    type: ["null", int]
    description: Number of samples that must pass the minPruning threshold
    inputBinding:
      position: 18
      prefix: "--numPruningSamples"

  - id: "#minReadsPerAlignmentStart"
    type: ["null", int]
    description: Minimum number of reads sharing the same alignment start for each genomic location in an active region
    inputBinding:
      position: 19
      prefix: "--minReadsPerAlignmentStart"

  - id: "#minPruning"
    type: ["null", int]
    description: Minimum support to not prune paths in the graph
    inputBinding:
      position: 20
      prefix: "--minPruning"

  - id: "#minDanglingBranchLength"
    type: ["null", int]
    description: Minimum length of a dangling branch to attempt recovery
    inputBinding:
      position: 21
      prefix: "--minDanglingBranchLength"

  - id: "#min_base_quality_score"
    type: ["null", int]
    description: Minimum base quality required to consider a base for calling
    inputBinding:
      position: 22
      prefix: "--min_base_quality_score"

  - id: "#maxReadsInRegionPerSample"
    type: ["null", int]
    description: Maximum reads in an active region
    inputBinding:
      position: 23
      prefix: "--maxReadsInRegionPerSample"

  - id: "#maxNumHaplotypesInPopulation"
    type: ["null", int]
    description: Maximum number of haplotypes to consider for your population
    inputBinding:
      position: 24
      prefix: "--maxNumHaplotypesInPopulation"
      
  - id: "#max_alternate_alleles"
    type: ["null", int]
    description: Maximum number of alternate alleles to genotype
    inputBinding:
      position: 25
      prefix: "--max_alternate_alleles"

  - id: "#kmerSize"
    description: Kmer size to use in the read threading assembler
    type:
      - "null"
      - type: array
        items: int
        inputBinding: { prefix: "--kmerSize " }
    inputBinding: 
      position: 26

  - id: "#input_prior"
    description: Input prior for calls
    type:
      - "null"
      - type: array
        items: double
        inputBinding: { prefix: "--input_prior" }
    inputBinding: 
      position: 27
    
  - id: "#ERCIS"
    type: ["null", int]
    description: The size of an indel to check for in the reference model
    inputBinding:
      position: 28
      prefix: "--indelSizeToEliminateInRefModel"

  - id: "#indel_heterozygosity"
    type: ["null", double]
    description: Heterozygosity for indel calling
    inputBinding:
      position: 29
      prefix: "--indel_heterozygosity"

  - id: "#heterozygosity"
    type: ["null", double]
    description: Heterozygosity for indel calling
    inputBinding:
      position: 30
      prefix: "--heterozygosity"

  - id: "#GVCFGQBands"
    description: Input prior for calls
    type:
      - "null"
      - type: array
        items: int
        inputBinding: { prefix: "--GVCFGQBands" }
    inputBinding: 
      position: 31

  - id: "#group"
    description: Input prior for calls
    type:
      - "null"
      - type: array
        items: string
        inputBinding: { prefix: "--group" }
    inputBinding: 
      position: 32

  - id: "#graphOutput"
    type: ["null", File]
    description: Write debug assembly graph information to this file
    inputBinding:
      position: 33
      prefix: "--graphOutput"

  - id: "#genotyping_mode"
    type: ["null", string]
    description: The --genotyping_mode argument is an enumerated type (GenotypingOutputMode), which can have one of the following values
    inputBinding:
      position: 34
      prefix: "--genotyping_mode"

  - id: "#gcpHMM"
    type: ["null", int]
    description: Flat gap continuation penalty for use in the Pair HMM
    inputBinding:
      position: 35
      prefix: "--gcpHMM"

  - id: "#forceActive"
    type: ["null", boolean]
    description: If provided, all bases will be tagged as active
    inputBinding:
      position: 36
      prefix: "--forceActive"

  - id: "#excludeAnnotation"
    description: One or more specific annotations to exclude
    type:
      - "null"
      - type: array
        items: string
        inputBinding: { prefix: "--excludeAnnotation" }
    inputBinding: 
      position: 37
      
  - id: "#emitRefConfidence"
    type: ["null", string]
    description: Mode for emitting reference confidence scores
    inputBinding:
      position: 38
      prefix: "--emitRefConfidence"

  - id: "#dontTrimActiveRegions"
    type: ["null", boolean]
    description: If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping
    inputBinding:
      position: 39
      prefix: "--dontTrimActiveRegions"

  - id: "#dontIncreaseKmerSizesForCycles"
    type: ["null", boolean]
    description: Disable iterating over kmer sizes when graph cycles are detected
    inputBinding:
      position: 39
      prefix: "--dontIncreaseKmerSizesForCycles"

  - id: "#doNotRunPhysicalPhasing"
    type: ["null", boolean]
    description: As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
    inputBinding:
      position: 40
      prefix: "--doNotRunPhysicalPhasing"

  - id: "#disableOptimizations"
    type: ["null", boolean]
    description: Dont skip calculations in ActiveRegions with no variants
    inputBinding:
      position: 41
      prefix: "--disableOptimizations"

  - id: "#debug"
    type: ["null", boolean]
    description: Print out very verbose debug information about each triggering active region
    inputBinding:
      position: 42
      prefix: "--debug"

  - id: "#contamination"
    type: ["null", File]
    description: Tab-separated File containing fraction of contamination in sequencing data (per sample) to aggressively remove. Format should be "" (Contamination is double) per line; No header.
    inputBinding:
      position: 43
      prefix: "--contamination_fraction_to_filter"
      
  - id: "#consensus"
    type: ["null", boolean]
    description: Print out very verbose debug information about each triggering active region
    inputBinding:
      position: 44
      prefix: "--consensus"

  - id: "#comp"
    description: comp binds reference ordered data. This argument supports ROD files of the following types BCF2, VCF, VCF3
    type:
      - "null"
      - type: array
        items: string
        inputBinding: { prefix: "--comp" }
    inputBinding: 
      position: 45

  - id: "#bandPassSigma"
    type: ["null", double]
    description: The sigma of the band pass filter Gaussian kernel; if not provided defaults to Walker annotated default
    inputBinding:
      position: 46
      prefix: "--consensus"

  - id: "#bamWriterType"
    type: ["null", string]
    description: Which haplotypes should be written to the BAM. 
    inputBinding:
      position: 47
      prefix: "--bamWriterType"
      
  - id: "#bamOutput"
    type: ["null", File]
    description: File to which assembled haplotypes should be written 
    inputBinding:
      position: 48
      prefix: "--bamOutput"

  - id: "#annotation"
    description: One or more specific annotations to apply to variant calls
    type:
      - "null"
      - type: array
        items: string
        inputBinding: { prefix: "--annotation" }
    inputBinding: 
      position: 49

  - id: "#annotateNDA"
    type: ["null", boolean]
    description: If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site 
    inputBinding:
      position: 50
      prefix: "--annotateNDA"

  - id: "#allSitePLs"
    type: ["null", boolean]
    description: Annotate all sites with PLs 
    inputBinding:
      position: 51
      prefix: "--allSitePLs"

  - id: "#allowNonUniqueKmersInRef"
    type: ["null", boolean]
    description: Allow graphs that have non-unique kmers in the reference 
    inputBinding:
      position: 52
      prefix: "--allowNonUniqueKmersInRef"

  - id: "#alleles"
    description: The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES
    type:
      - "null"
      - type: array
        items: string
        inputBinding: { prefix: "--alleles" }
    inputBinding: 
      position: 53

  - id: "#activityProfileOut"
    type: ["null", File]
    description: Output the raw activity profile results in IGV format 
    inputBinding:
      position: 54
      prefix: "--activityProfileOut"

  - id: "#activeRegionOut"
    type: ["null", File]
    description: Output the active region to this IGV formatted file 
    inputBinding:
      position: 55
      prefix: "--activeRegionOut"

  - id: "#activeRegionMaxSize"
    type: ["null", int]
    description: The active region maximum size; if not provided defaults to Walker annotated default 
    inputBinding:
      position: 56
      prefix: "--activeRegionMaxSize"

  - id: "#activeRegionExtension"
    type: ["null", int]
    description: The active region extension; if not provided defaults to Walker annotated default
    inputBinding: 
      position: 57
      prefix: "--activeRegionExtension"

  - id: "#activeProbabilityThreshold"
    type: ["null", double]
    description: Threshold for the probability of a profile state being active.
    inputBinding: 
      position: 58
      prefix: "--activeProbabilityThreshold"

outputs:

  - id: "#output_HaplotypeCaller"
    type: File
    outputBinding: 
      glob: $(inputs.outputfile_HaplotypeCaller)

arguments:
  - valueFrom: "./test/test-files"
    position: 2
    separate: false
    prefix: "-Djava.io.tmpdir="
    
  - valueFrom: "/usr/local/bin/GenomeAnalysisTK.jar"
    position: 3
    prefix: "-jar"

  - valueFrom: "HaplotypeCaller"
    position: 4
    prefix: "-T"

baseCommand: ["java"]


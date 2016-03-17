#!/usr/bin/env cwl-runner
cwlVersion: cwl:draft-3 
class: CommandLineTool
description: "Ray is a parallel de novo genome assembler that utilises the message-passing interface everywhere and is implemented using peer-to-peer communication."

$namespaces:
 edam: http://edamontology.org/
 foaf: http://xmlns.com/foaf/0.1/
 dct: http://purl.org/dc/terms/

$schemas:
 - http://edamontology.org/
 - http://dublincore.org/2012/06/14/dcterms.rdf
 - http://xmlns.com/foaf/spec/20140114.rdf

inputs:
 - id: kmer_length
   type: int
   description: "Selects the length of k-mers. The default value is 21."
   label: "k-mer length"
   inputBinding:
     prefix: "-k"
   
 - id: input_files
   label: "Paired end reads"
   description: "Provides two files containing paired-end reads."
   type:
     type: array
     items: File
     format: edam:format_1930
     inputBinding:
       prefix: "-p"
   
outputs:
 - id: contigs
   type: File
   description: "Contiguous sequences in FASTA format"
   format: edam:format_1929
   outputBinding:
    glob: "Contigs.fasta"
  
baseCommand: "Ray"

dct:creator:
- class: foaf:Organization
  foaf:name: "Norwegian University of Science and Technology"
  foaf:member:
  - class: foaf:Person
    foaf:mbox: "mailto:animesh.sharma@ntnu.no"
    foaf:name: "Animesh Sharma"
    
#/Ray -o testtt -p S2_1.fastq S2_2.fastq -k 31
#ContigLengths.txt
#Contigs.fasta
#Contigs.fasta.kmer.dna.output.default
#Contigs.fasta.kmer.dna.output.default.func.txt
#Contigs.fasta.kmer.dna.output.default.otu.txt
#CoverageDistributionAnalysis.txt
#CoverageDistribution.txt
#degreeDistribution.txt
#ElapsedTime.txt
#FilePartition.txt
#GraphPartition.txt
#LibraryData.xml
#LibraryStatistics.txt
#NeighbourhoodRelations.txt
#NetworkTest.txt
#NumberOfSequences.txt
#OutputNumbers.txt
#ParallelPaths.txt
#RayCommand.txt
#RayPlatform_Version.txt
#RaySmartCommand.txt
#RayVersion.txt
#ScaffoldComponents.txt
#ScaffoldLengths.txt
#ScaffoldLinks.txt
#Scaffolds.fasta
#SeedLengthDistribution.txt
#SequencePartition.txt

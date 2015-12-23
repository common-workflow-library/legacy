class: schema:SoftwareSourceCode
schema:name: "alea"
schema:about: >
  ALEA is a computational toolbox for allele-specific (AS) epigenomics analysis, which incorporates allelic variation data within existing
  resources, allowing for the identification of significant associations between epigenetic modifications and specific allelic variants
  in human and mouse cells. ALEA provides a customizable pipeline of command line tools for AS analysis of next-generation sequencing data
  (ChIP-seq, RNA-seq, etc.) that takes the raw sequencing data and produces separate allelic tracks ready to be viewed on genome browsers.
  ALEA takes advantage of the available genomic resources for human (The 1000 Genomes Project Consortium) and mouse (The Mouse Genome Project)
  to reconstruct diploid in silico genomes for human samples or hybrid mouse samples under study. Then, for each accompanying ChIP-seq or
  RNA-seq dataset, ALEA generates two wig files from short reads aligned differentially to each haplotype.
  This pipeline has been validated using human and hybrid mouse ChIP-seq and RNA-seq data (See Test Data section).

schema:targetProduct:
  class: schema:SoftwareApplication
  schema:downloadUrl: ftp://ftp.bcgsc.ca/supplementary/ALEA/files/alea.1.2.2.tar.gz
  schema:softwareVersion: "1.2.2"
  schema:applicationCategory: "command line tool"

schema:url: http://www.bcgsc.ca/platform/bioinfo/software/alea
schema:license: https://opensource.org/licenses/AFL-3.0
schema:programmingLanguage: "JAVA"

schema:publication:
- class: schema:ScholarlyArticle
  id: urn:pmid:24371156
  schema:name: "ALEA: a toolbox for allele-specific epigenomics analysis Bioinformatics (2014) 30 (8): 1172-1174."
  schema:author:
  - class: schema:Person
    schema:name: "Hamid Younesy"
  - class: schema:Person
    schema:name: "Torsten Moller"
  - class: schema:Person
    schema:name: "Alireza Heravi-Moussavi"
  - class: schema:Person
    schema:name: "Jeffrey B. Cheng"
  - class: schema:Person
    schema:name: "Joseph F. Costello"
  - class: schema:Person
    schema:name: "Matthew C. Lorincz"
  - class: schema:Person
    id: mailto:mkarimi@bcgsc.ca
  - class: schema:Person
    schema:name: "Steven J. M. Jones"
  schema:sameAs:
    - http://dx.doi.org/10.1093/bioinformatics/btt744
    - http://bioinformatics.oxfordjournals.org/content/30/8/1172.long

schema:author:
- class: schema:Person
  id: mailto:mkarimi@bcgsc.ca
  schema:name: "Mohammad Karimi"
  schema:email: mailto:mkarimi@bcgsc.ca
  schema:url: http://www.bcgsc.ca/author/mkarimi
  schema:worksFor:
  - class: schema:Organization
    schema:name: "Canada's Michael Smith Genome Sciences Centre, BC Cancer Agency, Vancouver, British Columbia, V5Z 4S6, Canada"
  - class: schema:Organization
    schema:name: "Department of Medical Genetics, Life Sciences Institute, The University of British Columbia, Vancouver, British Columbia, V6T 1Z3, Canada"


class: schema:SoftwareSourceCode
schema:name: "samtools"
schema:about: >
  A suite of programs for interacting with high-throughput sequencing data.
  It consists of three separate repositories: Samtools (Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format),
  BCFtools (Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants)
  and HTSlib (A C library for reading/writing high-throughput sequencing data).

schema:url: http://www.htslib.org/
schema:codeRepository: https://github.com/samtools/samtools.git

schema:license:
- https://opensource.org/licenses/MIT
- https://opensource.org/licenses/BSD-3-Clause

schema:targetProduct:
  class: schema:SoftwareApplication
  schema:softwareVersion: "1.2-242-4d56437"
  schema:applicationCategory: "command line tool"

schema:programmingLanguage: "C, Perl"

schema:publication:
- class: schema:ScholarlyArticle
  id: urn:pmid:21903627
  schema:name: "(Li, 2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics."
  schema:url: http://www.ncbi.nlm.nih.gov/pubmed/21903627
- class: schema:ScholarlyArticle
  id: urn:pmid:21320865
  schema:name: "(Li, 2011) Improving SNP discovery by base alignment quality. Bioinformatics."
  schema:url: http://www.ncbi.nlm.nih.gov/pubmed/21320865
- class: schema:ScholarlyArticle
  id: urn:pmid:19505943
  schema:name: "(Li et al., 2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics."
  schema:url: http://www.ncbi.nlm.nih.gov/pubmed/19505943

schema:discussionUrl:
- https://lists.sourceforge.net/lists/listinfo/samtools-help
- https://lists.sourceforge.net/lists/listinfo/samtools-devel

schema:creator:
- class: schema:Organization
  schema:name: "Sanger Institute"
  schema:member:
  - class: schema:Person
    schema:name: "Heng Li"
    schema:description: "wrote most of the initial source codes of SAMtools and various converters."
- class: schema:Organization
  schema:name: "Broad Institute"
  schema:member:
  - class: schema:Person
    schema:name: "Bob Handsaker"
    schema:description: |
      A major contributor to the
      SAM/BAM specification. He designed and implemented the BGZF format, the
      underlying indexable compression format for the BAM format. BGZF does
      not support arithmetic between file offsets.
- class: schema:Organization
  schema:name: "Beijing Genome Institute"
  schema:member:
  - class: schema:Person
    schema:name: "Jue Ruan"
    schema:description: |
      Designed and implemented the
      RAZF format, an alternative indexable compression format. RAZF is no longer
      used by or provided with SAMtools. Source code remains available in older
      SAMtools 0.1.x releases and from the standalone branch in the repository.
- class: schema:Person
  schema:name: "Colin Hercus"
  schema:description: "updated novo2sam.pl to support gapped alignment by novoalign."
- class: schema:Person
  schema:name: "Petr Danecek"
  schema:description: "contributed the header parsing library sam_header.c and sam2vcf.pl script."
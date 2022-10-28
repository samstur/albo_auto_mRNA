# Workflow/pipeline details for the albopictus autogeny RNA-seq analysis

## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found in the following google sheets: ([link](https://docs.google.com/spreadsheets/d/1_RixzDGNsUlvhOMVCTXViuR_Rim8dNrS5joLMOwbsXo/edit#gid=0))

### Data Accession
Data was generated at Georgetown by Mara, Peter, and Sam.

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control
Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/trim.sh))

FastQC (v0.11.9) was used for quality control visualization ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/fastqc.sh))

### Mapping
We will be mapping with STAR (v2.7.1a).

Mapping was done using the Aedes albopictus reference genome (aedes_albopictus_AalbF3.fa) and the corresponding annotation file (aedes_albopictus_AalbF3.gff3)

STAR (v2.7.1a) was used for indexing the genome ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/STAR_genomeIndex.sh))

Reads were mapped in a two pass method. The first pass followed typical method with splice junctions from annotations ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/STAR_mapping.sh)). The second pass is similar except that it additionally uses the output splice junctions info from the first pass (these would be novel splice junctions) to facilitate mapping ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/STAR_map_twopass.sh)).

Output sam files were converted to bam and then indexed ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/sam2bam.sh))

### HTSeq Gene Counts
HTSeq (v0.13.5) was used to counts reads mapped to genes for downstream analyses ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/htseq_count.sh))

# Downstream analysis of gene counts using R (4.1.3)

### DESeq package in R was used to identify differentially expressed genes
([script](https://github.com/samstur/albo_auto_mRNA/blob/main/DESeq_script.R))

### KEGGREST package in R was used to identify significantly enriched pathways
([script](https://github.com/samstur/albo_auto_mRNA/blob/main/KEGG_Enrichment.R))

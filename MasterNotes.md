# Workflow/pipeline details for the albopictus autogeny RNA-seq analysis

## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found in the following google sheets: ([link](https://docs.google.com/spreadsheets/d/1_RixzDGNsUlvhOMVCTXViuR_Rim8dNrS5joLMOwbsXo/edit#gid=0))

### Data Accession
Data was generated at Georgetown by Mara, Peter, and Sam.

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control
Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/trim.sh))

FastQC (v0.11.9) was used for quality control visualization ([script](https://github.com/samstur/albo_auto_mRNA/blob/main/fastqc.sh))

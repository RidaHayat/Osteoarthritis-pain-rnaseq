# Input data

The analysis uses the public bulk RNA-seq study of high- and low-pain knee osteoarthritis synovium: [NCBI BioProject PRJNA389258](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA389258) / [GEO GSE99662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99662). There are five samples in each pain group.

Download `GSE99662_gene_counts_matrix.txt.gz` from the GEO Series supplementary files and put it in this folder. The preparation script selects the 10 mRNA samples, maps their sample IDs and pain groups, and writes `counts_geo.tsv` and a minimal two-column `metadata.csv`. These local inputs are ignored by Git.

The supplied GEO file is a gene-level processed count matrix with fractional expected counts. No transcript-level quantifications or sample-specific effective lengths are supplied in this file, so a `tximport` workflow cannot be run from it. The analysis rounds fractional expected counts to the nearest integer for DESeq2. Do not substitute TPM or RPKM values for the gene-count matrix.


# Code
1. ExonArray: run_mmbgx.R: run MMBGX R package to generate the gene and isoform expressions from raw exon-array data.
2. Microarray: RMANomralization.m: matlab code to run RMA normalization to generate gene expressions from raw microarray gene expressiond data.
3. NanoString: (1) Counts_Normalized.mat: normalized raw NanoString data. (2) AMatrix: contains the NanoString indicator matrices (probe-by-isoform) for the genes. (3) GenerateIsoformExpression.m: matlab code to generate isoform expression from normalized raw NanoString data (equation 1 in the manuscript).
4. RNA-seq: RNA-seq_commands.txt: command lines to run both isoform and gene expressions quantification methods from raw RNA-seq fastq files.

# Data
1. NanoString: (1) IsoformExpression.xlsx: processed NanoString isoform expression data. (2) A.xlsx: ontains the NanoString indicator matrices (probe-by-isoform) for the genes. (3) 2ndReplicate.zip: a second technical replicate of raw NanoString data without house keeping genes.
2. RT-qPCR: RT-qPCR data for the four genes.

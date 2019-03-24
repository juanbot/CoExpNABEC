# CoExpNABEC

We have available in this package a co-expression network built from RNA-seq of 311 samples from frontal cortex, total RNA, stranded, aligned with STAR using hg19 UCSC reference genome and an initial quantification with salmon which quantifies at the transcript level with TPM (transcripts per million) quantification. 

Samples have been QCed and annotated in this regard. Also, there is exome sequencing and genotyping for a subset of these individuals. Most of these last subsets of individuals include methylation data.

The samples have been QCed for gender (based on methylation data and specific genes) and Ethnicity (based on population structure covariates from plink). Those samples appearing ats outliers were dropped off the project (there were only a few). 
The pipeline for RNA-seq and transcript quantification based on Salmon has been performed by Raph Gibbs at the NIH. Samples preparation has been coordinated by Mark Cookson also from the NIH.




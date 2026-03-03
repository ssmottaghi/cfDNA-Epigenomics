# cfDNA Analysis & cfDECOR Model Construction
This repository contains the end-to-end pipeline for cell-free DNA (cfDNA) analysis and the development of the cfDECOR model. The project is divided into two primary modules:

1. cfDNA-analysis: This folder contains the bioinformatics pipeline used to process cfDNA data and extract normalized depth at specified genomic regions.

methylation_analysis/: A specialized subfolder containing scripts and workflows to extract methylation values from Whole Genome Bisulfite Sequencing (WGBS) cfDNA data.


2. Model_Development: This folder contains the code for signature matrix construction and model training.

Signature Matrix: Processing of ATAC-Seq BED files to generate a reference signature matrix from chromatin accessibility data.

Deconvolution Model: The core development of the cfDECOR model, which utilizes the cfDNA depth extracted in the previous module for tissue-of-origin deconvolution.

README: Dataset for "Distinct microbiomes of the scleractinian coral Favia fragum in mangrove and adjacent reef habitats in the Panamanian Caribbean"

This repository contains the dataset and associated scripts for the manuscript "Distinct microbiomes of the scleractinian coral Favia fragum in mangrove and adjacent reef habitats in the Panamanian Caribbean." This dataset is associated with BioProject #PRJNA1023296 and SRA Submission #SUB13878240 of the National Center for Biotechnology Information.

Overview

This dataset comprises processed amplicon sequencing data from an Illumina MiSeq platform, along with a comprehensive workflow to reproduce the bioinformatic analyses presented in the associated manuscript. Users can access demultiplexed sequencing reads, follow a workflow for amplicon sequence variant (ASV) inference using DADA2, and perform downstream ecological analyses and visualizations.

Please note: Raw sequencing files are not provided in this Dryad repository. The demultiplexed Read 1 and Read 2 files per sample, which are the starting point for the DADA2 analysis in this study, are available directly from the Sequence Read Archive (SRA) under Submission #SUB13878240, for BioProject #PRJNA1023296.

Contents

This Dryad repository contains the following files:
	•	Favia_demultiplex_trim_workflow.txt: A generic and customizable workflow (text script) for demultiplexing and trimming raw Illumina MiSeq sequencing files. This script is provided for users to apply to their own raw sequencing data and to understand the initial processing steps that lead to demultiplexed reads.
	•	Favia_workflow_dada2_to_phyloseq.R: An R script to implement the DADA2 pipeline for amplicon analysis, generating a phyloseq object from demultiplexed reads.
	•	Favia.ps.RDS: A pre-computed phyloseq object in R's RDS format, representing the final processed amplicon dataset.
	•	Favia_workflow_phyloseq_analysis_visualization.R: An R script for analyzing community data, performing statistical tests, and generating visualizations from a phyloseq object.

Getting Started

Demultiplexing and Trimming Your Own Raw Data (Informational Only for this Dataset)

While the raw sequencing files for this specific study are not included here, we provide the Favia_demultiplex_trim_workflow.txt script as a generic guide for demultiplexing and trimming Illumina MiSeq raw data. This workflow outlines the steps that produce demultiplexed Read 1 and Read 2 files per sample, consistent with the format of the files available on SRA for this study. You can use this script to process your own raw sequencing data or to understand the initial steps taken in similar amplicon sequencing pipelines.
	1	Download Raw Data: Obtain your raw Illumina MiSeq sequencing files (e.g., FASTQ files).
	2	Demultiplexing and Trimming Workflow:
		◦	Open and customize the Favia_demultiplex_trim_workflow.txt script.
		◦	This workflow requires the Perl script "MergeMeCheck3.pl". This script is available on our GitHub website: https://marineinvert.github.io/microbiome/resources.html. Please download MergeMeCheck3.pl and place it in an accessible directory.
		◦	Execute the steps outlined in Favia_demultiplex_trim_workflow.txt using your preferred command-line environment.
	3	Outputs: The output of this workflow will be demultiplexed Read 1 and Read 2 FASTQ files for each sample, along with statistics from the demultiplexing and trimming steps.


Users have two primary options for working with this dataset:
	1	Start from demultiplexed files (downloaded from SRA): Use the demultiplexed Read 1 and Read 2 files from SRA as the input for the DADA2 pipeline.
	2	Start from the supplied phyloseq object: Directly use the pre-computed Favia.ps.RDS phyloseq object for downstream community analyses.



1.1 Amplicon Analysis with DADA2 (from Demultiplexed Files)

The demultiplexed Read 1 and Read 2 files for this study are available in the Sequence Read Archive (SRA) under Submission #SUB13878240 for BioProject #PRJNA1023296. If you have downloaded these files, you can proceed with ASV inference using DADA2:
	1	Download Demultiplexed Files: Obtain the demultiplexed Read 1 and Read 2 FASTQ files for each sample from SRA (BioProject #PRJNA1023296, Submission #SUB13878240).
	2	Download and Install R and RStudio: Ensure you have a recent version of R and RStudio installed.
	3	Install Required R Packages: Open RStudio and install the necessary packages, particularly dada2 and phyloseq.
	4	Run Favia_workflow_dada2_to_phyloseq.R:
		◦	Open the Favia_workflow_dada2_to_phyloseq.R script in RStudio.
		◦	Download the Favia_metadata_2.csv table containing the metadata for each sample and save it into an accessible directory.
		◦	Modify the file paths within the script to point to your downloaded demultiplexed Read 1 and Read 2 FASTQ files.
		◦	Run the script. This script implements the DADA2 pipeline, including quality filtering, dereplication, ASV inference, chimera removal, and taxonomic assignment.
		◦	Outputs of DADA2: The script will generate a phyloseq object. This object encapsulates:
			▪	ASV Table: A matrix of amplicon sequence variants (ASVs) by samples, containing the abundance of each ASV in each sample.
			▪	Taxonomy Table: A table assigning taxonomic classifications (e.g., Kingdom, Phylum, Class, Order, Family, Genus, Species) to each ASV.
			▪	Sample Data: A data frame containing metadata for each sample (e.g., experimental conditions, sampling location, etc.).
			▪	Phylogenetic Tree (Optional): Depending on the specific DADA2 pipeline implementation, a phylogenetic tree relating the ASVs may also be included.

1.2 Community Analysis and Visualization (from Phyloseq Object)

You can either use the phyloseq object generated from the previous step or directly use the supplied Favia.ps.RDS file to perform downstream analyses:
	1	Load Phyloseq Object:
		◦	If you generated your own phyloseq object, it will be in your R environment.
		◦	If using the supplied object, load it into R using readRDS("Favia.ps.RDS").
	2	Run Favia_workflow_phyloseq_analysis_visualization.R:
		◦	Open the Favia_workflow_phyloseq_analysis_visualization.R script in RStudio.
		◦	Run the script to reproduce the analyses and visualizations.

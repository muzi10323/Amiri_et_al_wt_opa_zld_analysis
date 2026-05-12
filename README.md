This repository contains the R-based code for ATAC-seq and RNAn-seq data analysis in Aimir et al 2026. 

📁 Repository Structure
Rmd/: Contains the core R scripts and R Markdown files for data importation, processing, and visualization.

Import_Clogmia_Paired_End_bam_RNAseq_to_R_v1.R: Scripts for processing and importing paired-end RNA-seq data.

Import_Clogmia_Paired_End_bam_to_R_v3.R: Scripts for processing and importing paired-end ATAC-seq data.

26_ALL_Analysis_8_clean.Rmd: The primary, finalized analysis script integrating the datasets.

250401_New_Wiggles.Rmd: Script for generating wiggle tracks/coverage plots for genome browser visualization.

26_ALL_Analysis_8_clean_files/: Automatically generated figures and cached output from the main Rmd analysis.

PWMs/: Position Weight Matrices used for transcription factor binding motif analysis.

Experimental Conditions:

WT_timecourse/: Output and specific data subsets (including count.txt files) for the wild-type developmental timecourse.

zldKD/: Analysis files and counts specific to the Zelda knockdown experiments.

opaKD/: Analysis files and counts specific to the Opa knockdown experiments.

Root Files:

all-coldata.xlsx: The master metadata file detailing conditions, batches, and sample information for the integrated epigenomic and transcriptomic datasets.

BSgenome.Calbipunctata...: Custom compiled BSgenome package required for sequence extraction and motif mapping in the Clogmia albipunctata genome.

26_ALL_Analysis_8_clean.html: The compiled HTML report of the finalized analysis.

🚀 Usage & Reproduction
1. Environment Setup
Genome Package: Before running any scripts, you must install the custom BSgenome.Calbipunctata package. Install it via terminal using R CMD INSTALL or within R using:

R
install.packages("BSgenome.Calbipunctata...", repos = NULL, type = "source")
2. Data Import Options
You can reproduce the analysis either from raw sequencing files or by using the provided pre-computed count files.

Option A: From Raw Sequencing Files

Process the raw reads using the pipeline outlined in the methods section of the manuscript.

Import the processed raw read files using the Import_Clogmia_Paired_End_... R scripts located in the Rmd/ folder. Ensure you map the files according to the metadata in all-coldata.xlsx.

Store the newly imported ATAC-seq and RNA-seq data in a directory named ./GRanges to replace the empty sample files.

Option B: From Pre-computed Counts
Alternatively, you can skip the raw read importation by using the count.txt files provided in the WT_timecourse/, opaKD/, and zldKD/ folders for each respective experiment.

Note: If you choose this route, you will need to modify the Rmd code to read directly from these text files rather than looking for the imported GRanges objects.

3. Execution
Run the main analysis scripts located in the Rmd/ directory. It is highly recommended to "knit" 26_ALL_Analysis_8_clean.Rmd—this will seamlessly execute the code, reproduce the full HTML report, and automatically update the output figures in the 26_ALL_Analysis_8_clean_files/ directory.

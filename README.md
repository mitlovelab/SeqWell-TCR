### Tu, A. A. et al. TCR sequencing paired with massively parallel 3′ RNA-seq reveals clonotypic T cell signatures. Nat. Immunol. 2019 2012 20, 1692–1699 (2019).
Landing page for computational scripts and other supplementary information related to TCR recovery from Seq-Well 
libraries related to the Nat. Immunol manuscript.

#### GEO
Gene expression matrices can be found on GSE136028. Seurat objest will be included as rds files.

#### Rscripts

Contains basic R scripts used to produce Figure 3 and 4 in the manuscript. Also includes tests and other exploratory
analysis that ultimately did not get included.

##### Figure 3
The E7_combined script contains analysis that compares TCRs from mice that were immunized with E7, and sorted with 
tetramers. 

##### Figure 4
PNOIT_TCR script contains analysis regarding peanut allergy patients shown in figure 4.

##### src
Contains supplemental files (such as gene lists, scripts for helper functions), as well as msigdb results that were 
obtained via web portal. The `MsigDb` folder contains results obtained from querying MsigDb portal. Include data for both
figure 3 and 4.

##### Preprocess
Contains script for combining raw TCR data (as made available in supplement of manuscript) with the relevant seurat files

#### TCR_fastqprocessing
Shell scripts for mapping sequencing fastq to TCR reference to produce the relevant bam files. Reference files included
for both human and mouse

#### Seqwell_TCR_processing
Instruction for analyzing mapped bam files into compiled TCR summary files, which contains TCR sequences along with the
single-cell barcodes that they are mapped to. Implemented as a docker image. Raw scripts included in this folder

##### Local_processing
Running the processing program locally just using MatLab. This is not recommended, as it can take a really long time.

##### Docker_associated_scripts
Scripts included and used in the docker image environment

#### TCR_recovery_protocol
Laboratory protocol on how to perform the recovery procedure to construct the TCR sequencing libraries.

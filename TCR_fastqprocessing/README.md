## Fastq Processing of TCR sequencing data
Scripts related to processing raw sequencing data from sequencing the TCR-enriched libraries. Currently, the code runs on Linux Server Luria provided by the MIT BioMicroCenter, which uses slurm job manager.
The Scripts will likely require other modifications to run on a different distribution of the linux server.

### Overview
The scripts here pre-process the raw sequencing fastqs, which is done with single-ended Illumina sequencing. Read1 targets the CDR3 region, and the index read targets the cell barcodes. 

The scripts (.sh files) takes sequencing data,and first rearrange the sequences into a standard read1 (cell barcode) and read2 (CDR3 mapping reads) format. The scripts then map it to the generated TCR reference files (provided in the respetive folders) in mouse and human via BWA, resulting in a single Bam file (TCRalignSort). The Bam file is then used to pass into analysis (most commonly run on cloud instance via docker). Separate human and mouse scripts are provided.

### Notes
This step was originally used to filter reads that do not map to any of the TCR regions. With the final enrichment protocol, not many reads are typically filtered out. This step may be removed in future updates.


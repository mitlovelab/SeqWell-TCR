### Seq-Well TCR processing and analysis

After fastq processing, the resulting bam file is used to make clonotype call on a per cell barcode/UMI basis. This page is the landing
page for TCR processing code and Docker documentation.

#### Docker
The easiest way to run the analysis is through Docker, ideally on a cloud (we use google cloud) instance. Please see 
https://hub.docker.com/orgs/mitlovelab/repositories

Also see the pdf instruction.

We are currently working on migrating the workflow onto Broad Terra, and into a more efficient combination of WDL and Python. Stay tune to this github for updates. 

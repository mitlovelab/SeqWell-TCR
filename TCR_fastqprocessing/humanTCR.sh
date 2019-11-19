#!/bin/sh
# request Bourne shell as shell for job
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m be
#$ -pe whole_nodes 1
###################  

#The goal of the script is to first transform the Read-1-index-1 format to Read-1-Read-2 format. that way, the TCR CDR3
#sequence should be in Read 1, and the associated cell barcode/UMI should be in Read 2.
#Then, the resulting pair of Fastqs are mapped to TCR reference (via bwa) to filter for reads that map (anywhere) on the TCR
#loci (the referecne is a collection of all VDJ and C genes). The final output is a bam file, which will then be used for subsequent processing. 
#Please see submission syntax example below for how to kick off this script in via bash.

module load bwa/0.7.12
module load samtools/1.3

awk 'NR%4==2' "$1" | rev | tr ATCG TAGC > "$2_IndexReads.txt"
awk 'NR%4==0' "$1" > "$2_QSeq.txt"
awk 'NR%4==1' "$1" | grep -o "[ATCGN]*" > "$2_BCSeq.txt"
sed 's~[ATGCN]~@~g' "$2_BCSeq.txt" > "$2_QRead2.txt"
awk 'NR%4==1' "$1" | grep -o "^.*#" > "$2_seqHeaders.txt" 
sed 's~#~#/1~' "$2_seqHeaders.txt" > "$2_seqHeadersRead1.txt"
sed 's~#~#/2~' "$2_seqHeaders.txt" > "$2_seqHeadersRead2.txt"
sed 's~@~+~' "$2_seqHeadersRead2.txt" > "$2_qualHeadersRead2.txt"
sed 's~@~+~' "$2_seqHeadersRead1.txt" > "$2_qualHeadersRead1.txt"
paste -d '\n' "$2_seqHeadersRead2.txt" "$2_IndexReads.txt" "$2_qualHeadersRead2.txt" "$2_QSeq.txt" > "$2_TCRreadFinal.fastq"
paste -d '\n' "$2_seqHeadersRead1.txt" "$2_BCSeq.txt" "$2_qualHeadersRead1.txt" "$2_QRead2.txt" > "$2_Read1.fastq"

#first variable is the sequence we put in. Second one is the sample name. third is the user name
#example: qsub humanTCR.sh fastq outputPrefix user

bwa mem /home/$3/data/humanTCR/humanTCR.fa $2_Read1.fastq $2_TCRreadFinal.fastq | samtools view -hF 256 > $2_TCRalign.sam
samtools sort -o $2_TCRalignSort.bam $2_TCRalign.sam
samtools index $2_TCRalignSort.bam

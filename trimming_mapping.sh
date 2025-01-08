#!/bin/bash

PREFIX=$1
module load conda
module load trimmomatic
module load bwa
module load samtools

# Define input files based on the provided prefix
READ1=${PREFIX}_1.fq.gz
READ2=${PREFIX}_2.fq.gz

# Create output directories if they don't exist
mkdir -p 1_fastq bam

# Step 1: Trim the reads using Trimmomatic
java -jar $TRIMMOMATIC/trimmomatic.jar PE -threads 8 -phred33 \
$READ1 $READ2 \
1_fastq/${PREFIX}_1.trimmed.fastq.gz 1_fastq/${PREFIX}_1un.trimmed.fastq.gz \
1_fastq/${PREFIX}_2.trimmed.fastq.gz 1_fastq/${PREFIX}_2un.trimmed.fastq.gz \
SLIDINGWINDOW:4:20

# Define the trimmed reads
TRIMMED_READ1=1_fastq/${PREFIX}_1.trimmed.fastq.gz
TRIMMED_READ2=1_fastq/${PREFIX}_2.trimmed.fastq.gz

# Define the reference genome
REF=10.fasta

# Step 2: Align the reads to the reference genome using BWA
bwa mem -t 32 -R "@RG\tID:$PREFIX\tSM:$PREFIX\tPL:ILLUMINA" $REF $TRIMMED_READ1 $TRIMMED_READ2 | \
samtools sort -n -@5 -o bam/${PREFIX}.bam

# Step 3: Fix mate information
samtools fixmate -m bam/${PREFIX}.bam - | \
samtools sort -@5 -o bam/${PREFIX}.fix.bam

# Step 4: Index the fixed BAM file
samtools index bam/${PREFIX}.fix.bam

# Step 5: Mark duplicates
samtools markdup -@5 -s bam/${PREFIX}.fix.bam - | \
samtools sort -@5 -o bam/${PREFIX}.fix.markdup.bam

# Step 6: Index the marked duplicates BAM file
samtools index bam/${PREFIX}.fix.markdup.bam

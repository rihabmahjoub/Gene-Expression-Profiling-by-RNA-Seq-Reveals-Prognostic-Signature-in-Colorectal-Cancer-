#!/bin/bash

# RNA-seq Pipeline - Part 1: Bash Shell Script
# Author: Rihab Mahjoub
# Description: Downloads RNA-seq data, performs quality control, filtering, and quantification using kallisto

# 1. Create folders
mkdir -p fastq_files filtered_fastq fastqc_dir kallisto_output Ref
cd fastq_files

# 2. Download FASTQ files manually if necessary (example for one sample)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR975/001/SRR975551/SRR975551_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR975/001/SRR975551/SRR975551_2.fastq.gz
# Repeat for SRR975552, SRR975569, SRR975570 if not already done

# 3. Go to Ref directory and download the transcriptome index
cd ../Ref
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# 4. Build kallisto index
kallisto index -i HS_Kallisto_index Homo_sapiens.GRCh38.cdna.all.fa.gz

# 5. Perform quality control with FastQC
cd ../fastq_files
fastqc *.fastq.gz -o ../fastqc_dir -t 4

# 6. Run MultiQC on FastQC results
cd ../fastqc_dir
multiqc .

# 7. Filter reads with fastp
cd ../fastq_files
for SAMPLE in SRR975551 SRR975552 SRR975569 SRR975570
  do
    fastp \
    -i ${SAMPLE}_1.fastq.gz \
    -I ${SAMPLE}_2.fastq.gz \
    -o ../filtered_fastq/${SAMPLE}_1.filtered.fastq.gz \
    -O ../filtered_fastq/${SAMPLE}_2.filtered.fastq.gz \
    -q 30
  done

# 8. Run Kallisto quantification
cd ../filtered_fastq
for SAMPLE in SRR975551 SRR975552 SRR975569 SRR975570
  do
    kallisto quant \
    -i ../Ref/HS_Kallisto_index \
    -o ../kallisto_output/${SAMPLE} \
    ${SAMPLE}_1.filtered.fastq.gz ${SAMPLE}_2.filtered.fastq.gz
  done

# RNA-Seq Pipeline: Tumor vs Normal (Kallisto + DESeq2 + Enrichment)

## Overview
This project performs RNA-seq quantification using **Kallisto** followed by differential expression and enrichment analysis using **R**.

### 📁 Structure
- `code.sh`: Bash script for FASTQ download, filtering, and quantification.
- `code2.r`: R script for differential expression, GO, and GSEA.
- `kallisto_output/`: Contains Kallisto quantification results.
- `fastqc_dir/`, `filtered_fastq/`, etc.: Quality control and filtered data.
- `Ref/`: Transcriptome reference used for Kallisto indexing.



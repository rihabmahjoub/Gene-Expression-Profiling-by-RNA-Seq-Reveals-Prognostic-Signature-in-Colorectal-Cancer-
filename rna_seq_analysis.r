# RNA-seq Analysis - Part 2: R Script
# Author: Rihab Mahjoub
# Description: Gene-level summarization, differential expression, GO and GSEA enrichment

# ---- Load Required Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("tximport", "EnsDb.Hsapiens.v86", "DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))
install.packages(c("readr", "dplyr", "tibble", "ggplot2", "cowplot", "matrixStats"))

library(tximport)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(matrixStats)

# ---- Load Kallisto Quantifications ----
samples_df <- read_tsv("kallisto_output/study_design.txt")
files <- file.path("kallisto_output", samples_df$sample, "abundance.tsv")
names(files) <- samples_df$sample

Tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name")) %>%
  as_tibble() %>%
  rename(target_id = tx_id) %>%
  select(target_id, gene_name)

txi_gene <- tximport(files,
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM")

# ---- Differential Expression with DESeq2 ----
countData <- round(txi_gene$counts)
sample <- data.frame(
  row.names = colnames(countData),
  group = c("tumor", "tumor", "normal", "normal")
)

dds <- DESeqDataSetFromMatrix(countData, sample, design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "tumor", "normal"))

# ---- Volcano Plot ----
res_df <- as.data.frame(na.omit(res))
res_df$gene <- rownames(res_df)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text(data = subset(res_df, significant == "Significant"), aes(label = gene), size = 2, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Padj")

# ---- GO Enrichment (gprofiler2) ----
BiocManager::install("gprofiler2")
library(gprofiler2)
resSig <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
gost_res <- gost(rownames(resSig), organism = "hsapiens")
gostplot(gost_res)

# ---- GSEA ----
C3 <- read.gmt("c4.all.v2025.1.Hs.symbols.gmt")
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$gene
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- GSEA(gene_list, TERM2GENE = C3, pvalueCutoff = 0.05)

# Plot
gseaplot(gsea_res, geneSetID = 1)
dotplot(gsea_res, showCategory = 10)

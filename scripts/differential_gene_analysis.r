# =============================================================================
# Automated RNA-seq Differential Expression Analysis Pipeline
# =============================================================================
#
# Project:
# Comparative Transcriptomic Analysis of Experimental Conditions
#
# Description:
# This R workflow performs automated differential gene expression analysis
# using DESeq2 on gene-level count data generated from RNA-seq experiments.
# The pipeline imports featureCounts output, constructs a count matrix,
# performs normalization and statistical testing, and generates publication-
# ready visualizations for exploratory and differential expression analysis.
#
# The workflow identifies significantly differentially expressed genes
# between experimental conditions across multiple time points and produces
# quality control and downstream analysis plots for biological interpretation.
#
# Workflow:
#   1. Import featureCounts gene count matrix
#   2. Create sample metadata table
#   3. DESeq2 normalization and dispersion estimation
#   4. Differential expression analysis
#   5. Variance stabilizing transformation (VST)
#   6. Principal Component Analysis (PCA)
#   7. Volcano plot generation
#   8. MA plot generation
#   9. Heatmap visualization of significant genes
#
# Input:
#   - gene_counts.txt (featureCounts output)
#
# Output:
#   - DEG_Day4.csv
#   - DEG_Day6.csv
#   - PCA_plot.pdf
#   - Volcano_Day4.pdf
#   - Volcano_Day6.pdf
#   - MA_Day4.pdf
#   - MA_Day6.pdf
#   - Heatmap_Top50_DEGs.pdf
#
# R Packages:
#   - DESeq2
#   - EnhancedVolcano
#   - pheatmap
#   - apeglm
#
# Author: Anshu Mathuria
# MSc Biotechnology & Bioinformatics
# Institute of Bioinformatics and Applied Biotechnology (IBAB)
#
# =============================================================================
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(apeglm)

# =====================================================
# READ FEATURECOUNTS OUTPUT
# =====================================================

counts <- read.delim(
  "gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

gene_ids <- counts$Geneid

count_matrix <- counts[,7:ncol(counts)]

rownames(count_matrix) <- gene_ids

# =====================================================
# SAMPLE METADATA
# MODIFY SAMPLE NAMES TO MATCH YOUR DATA
# =====================================================

sample_info <- data.frame(
  row.names = colnames(count_matrix),

  Condition = c(
    "WT_D4","WT_D4","WT_D4",
    "MUT_D4","MUT_D4","MUT_D4",
    "WT_D6","WT_D6","WT_D6",
    "MUT_D6","MUT_D6","MUT_D6"
  )
)

# =====================================================
# CREATE DESEQ OBJECT
# =====================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ Condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

# =====================================================
# VARIANCE STABILIZATION
# =====================================================

vsd <- vst(dds, blind = FALSE)


# =====================================================
# DAY 4 COMPARISON
# WT_D4 vs MUT_D4
# =====================================================

res_d4 <- results(
  dds,
  contrast = c(
    "Condition",
    "MUT_D4",
    "WT_D4"
  )
)

write.csv(
  as.data.frame(res_d4),
  "DEG_Day4.csv"
)

pdf(
  "Volcano_Day4.pdf",
  width = 8,
  height = 6
)

EnhancedVolcano(
  res_d4,
  lab = rownames(res_d4),
  x = "log2FoldChange",
  y = "padj",
  title = "WT Day4 vs MUT Day4"
)

dev.off()

pdf(
  "MA_Day4.pdf",
  width = 8,
  height = 6
)

plotMA(
  res_d4,
  ylim = c(-5,5)
)

dev.off()

# =====================================================
# DAY 6 COMPARISON
# WT_D6 vs MUT_D6
# =====================================================

res_d6 <- results(
  dds,
  contrast = c(
    "Condition",
    "MUT_D6",
    "WT_D6"
  )
)

write.csv(
  as.data.frame(res_d6),
  "DEG_Day6.csv"
)

pdf(
  "Volcano_Day6.pdf",
  width = 8,
  height = 6
)

EnhancedVolcano(
  res_d6,
  lab = rownames(res_d6),
  x = "log2FoldChange",
  y = "padj",
  title = "WT Day6 vs MUT Day6"
)

dev.off()

pdf(
  "MA_Day6.pdf",
  width = 8,
  height = 6
)

plotMA(
  res_d6,
  ylim = c(-5,5)
)

dev.off()

# =====================================================
# TOP 50 GENES HEATMAP (DAY 6)
# =====================================================

sig_genes <- rownames(
  subset(
    as.data.frame(res_d6),
    padj < 0.05
  )
)

if(length(sig_genes) > 50){
  sig_genes <- sig_genes[1:50]
}

heatmap_matrix <- assay(vsd)[sig_genes,]

pdf(
  "Heatmap_Top50_DEGs.pdf",
  width = 8,
  height = 10
)

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = sample_info,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_rownames = TRUE
)

dev.off()

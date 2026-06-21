# Methodology

## Project Overview

This project combines transcriptomic analysis, viral genome mapping, chimeric transcript detection, differential gene expression analysis, and experimental validation to investigate host–virus interactions in *Saccharomyces cerevisiae*.

---

## Part I: Host–Virus Chimeric Transcript Detection Workflow

### 1. RNA-Seq Quality Assessment

Raw paired-end RNA-seq datasets from 48 samples were assessed using FastQC. Quality metrics including read quality scores, GC content, sequence duplication levels, and adapter contamination were evaluated before downstream analysis.

### 2. Host Genome Alignment

High-quality reads were aligned to the *Saccharomyces cerevisiae* reference genome using STAR. Alignment statistics were collected to evaluate mapping efficiency across all samples.

### 3. Extraction of Unmapped Reads

Reads that failed to align to the host genome were extracted and retained for viral genome analysis.

### 4. Viral Genome Mapping

Unmapped reads were aligned against reference viral genomes to identify virus-associated transcripts and estimate viral abundance across samples.

### 5. Chimeric Alignment

Combined host–virus references were used for chimeric alignment using STAR's chimeric mapping functionality to identify candidate host–viral junctions.

### 6. Junction Analysis

STAR-generated chimeric junction files were processed using custom Python workflows. Candidate junctions were filtered based on breakpoint information and duplicate events were removed to improve confidence.

### 7. Fusion Sequence Reconstruction

BioPython-based scripts were used to reconstruct candidate fusion sequences by extracting genomic regions surrounding host and viral breakpoints.

### 8. Data Visualization

Python-based analyses were used to generate mapping summaries, read distribution plots, junction statistics, and exploratory visualizations.

### 9. Experimental Validation

Fusion sequences were used for primer design using Primer3Plus. Selected candidates were validated through RT-qPCR experiments.

---

## Part II: Differential Expression Analysis

### 1. Read Count Generation

Gene-level read counts generated from host genome alignments were used as input for transcriptomic analysis.

### 2. Differential Expression Analysis

Differential expression analysis was performed in R using DESeq2. Multiple biological comparisons were analyzed to identify genes showing significant expression changes.

### 3. Statistical Analysis

Genes were evaluated using normalized counts, log2 fold changes, adjusted p-values, and false discovery rate (FDR) correction.

### 4. Functional Interpretation

Differentially expressed genes were examined to identify transcriptomic signatures associated with host cellular responses.

### 5. Visualization

Results were visualized using:

* Volcano Plots
* Heatmaps
* Expression Profiles
* Differential Expression Summaries

---

## Computational Skills Demonstrated

* RNA-Seq Analysis
* Transcriptomics
* Differential Gene Expression Analysis
* Viral Genomics
* Chimeric Read Detection
* Genome Alignment
* Python-Based Data Analysis
* R-Based Statistical Analysis
* Data Visualization
* Bioinformatics Pipeline Development
* Linux Workflow Automation
* Experimental Validation Design

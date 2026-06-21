# Pipeline Logic

## Overview

The objective of this workflow was to identify potential host–viral chimeric transcripts and investigate transcriptomic changes associated with viral persistence using RNA sequencing data.

The pipeline was designed to progressively reduce data complexity while retaining biologically relevant information for downstream analysis and experimental validation.

---

## Step 1: Quality Assessment

Before any analysis, raw sequencing reads were assessed to ensure sufficient quality for reliable alignment and downstream interpretation.

Poor-quality reads, sequencing artifacts, or adapter contamination can introduce mapping errors and false-positive biological signals.

**Purpose:** Establish data quality and identify potential technical issues.

---

## Step 2: Host Genome Alignment

Reads were first aligned to the Saccharomyces cerevisiae reference genome.

Since the majority of RNA-seq reads originate from the host transcriptome, this step allows accurate quantification of host gene expression while separating reads that cannot be explained by the host genome.

**Purpose:** Characterize host transcriptional activity and identify unmapped reads for further investigation.

---

## Step 3: Unmapped Read Extraction

Reads that failed to align to the host genome were isolated.

These reads may represent viral transcripts, novel sequences, low-complexity regions, sequencing artifacts, or potential host–viral junctions.

Rather than discarding them, they were retained for targeted analysis.

**Purpose:** Enrich for candidate non-host sequences.

---

## Step 4: Viral Genome Mapping

Unmapped reads were aligned against viral reference genomes.

This step determines whether a proportion of the unexplained reads originate from viral RNA and provides an estimate of viral transcript abundance across samples.

**Purpose:** Identify virus-associated reads and characterize viral activity.

---

## Step 5: Chimeric Alignment

To investigate potential interactions between host and viral transcripts, combined host–virus references were used for chimeric alignment.

STAR chimeric mapping reports reads whose segments align to different genomic locations or reference sequences.

These events can indicate candidate host–viral junctions.

**Purpose:** Detect potential host–viral chimeric transcripts.

---

## Step 6: Junction Processing and Filtering

Raw chimeric junction outputs may contain alignment artifacts, low-confidence events, or duplicate observations.

Custom Python workflows were used to process junction files, remove redundant events, and prioritize high-confidence candidates.

**Purpose:** Improve confidence and reduce false-positive detections.

---

## Step 7: Fusion Sequence Reconstruction

For each candidate event, genomic coordinates were used to reconstruct the predicted fusion sequence.

BioPython-based approaches enabled extraction and assembly of sequence regions surrounding junction breakpoints.

These reconstructed sequences served as templates for downstream validation.

**Purpose:** Generate biologically interpretable fusion candidates.

---

## Step 8: Visualization and Data Exploration

Mapping statistics and candidate junction summaries were visualized using Python-based analytical workflows.

Visualization facilitates quality assessment, pattern discovery, and communication of findings.

**Purpose:** Interpret complex sequencing datasets and identify trends.

---

## Step 9: Experimental Validation

Candidate fusion sequences were used for primer design using Primer3Plus.

Selected events were subjected to RT-qPCR validation to evaluate whether computationally predicted candidates could be experimentally detected.

**Purpose:** Provide independent validation of computational observations.

---

# Differential Expression Analysis Logic

The host-aligned RNA-seq data was simultaneously used for transcriptomic analysis.

---

## Read Count Analysis

Gene-level count data generated from host genome alignments were used as input for differential expression analysis.

This enables quantitative comparison of transcriptional activity between biological conditions.

**Purpose:** Measure changes in host gene expression.

---

## DESeq2 Statistical Framework

Differential expression analysis was performed using DESeq2 in R.

Normalization, dispersion estimation, and statistical testing were applied to identify significantly altered genes.

**Purpose:** Detect biologically meaningful transcriptional changes.

---

## Visualization and Interpretation

Differential expression results were visualized using:

* Volcano plots
* Heatmaps
* Gene expression summaries

These visualizations facilitate identification of significantly regulated genes and broader transcriptomic trends.

**Purpose:** Interpret host cellular responses and generate biological hypotheses.

---

# Overall Workflow Strategy

The pipeline follows two complementary analytical branches:

### Branch 1: Host–Virus Interaction Analysis

RNA-seq → Host Alignment → Viral Mapping → Chimeric Detection → Fusion Reconstruction → RT-qPCR Validation

### Branch 2: Host Transcriptome Analysis

RNA-seq → Host Alignment → Read Counts → DESeq2 → Differential Expression Analysis → Biological Interpretation

Together, these approaches provide both a targeted investigation of potential host–viral interactions and a global assessment of host transcriptional responses.


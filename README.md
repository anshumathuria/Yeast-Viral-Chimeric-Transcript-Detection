# RNA-Seq Host–Virus Chimeric Transcript Detection Pipeline

## Overview

This repository presents the computational workflow developed as part of a Master's thesis project investigating host–virus interactions and transcriptomic signatures associated with persistent viral replication in *Saccharomyces cerevisiae* using RNA-Seq data.

The study focused on identifying and characterizing potential host–viral chimeric transcripts through sequential genome alignment, viral mapping, chimeric junction detection, and downstream bioinformatics analyses. The repository serves as a demonstration of the computational methodology, tools, and analytical approaches employed in the project.

## Research Focus

Understanding the interplay between viral persistence, transcriptome dynamics, and genomic stability through large-scale RNA sequencing and computational analysis.

## Dataset

* Organism: *Saccharomyces cerevisiae*
* High-throughput paired-end RNA-Seq datasets
* Multiple biological conditions and time-course samples
* Total Samples Analysed: 48

## Computational Workflow

1. Raw RNA-Seq Quality Assessment
2. Reference Genome Alignment
3. Differential expression of genes
4. Extraction of Unmapped Reads
5. Viral Genome Mapping
6. Combined Host–Virus Reference Construction
7. Chimeric Read Detection
8. Junction Generation and Filtering
9. Candidate Fusion Analysis
10. Validation-Oriented Mapping and Quality Assessment
11. Data Visualization and Interpretation
12. RT-qPCR Validation of the Chimeric reads

## Tools and Software

### RNA-Seq Processing

* FastQC
* STAR
* Bowtie2
* SAMtools

### Data Analysis and Statistics

* Python
* Pandas
* NumPy
* R

### Visualization

* Matplotlib
* Seaborn
* ggplot2

### Computing Environment

* Linux
* Bash
* Conda

## Skills Demonstrated

* RNA-Seq Data Analysis
* Transcriptomics
* Genome Alignment and Mapping
* Viral Genomics
* Chimeric Read Detection
* Bioinformatics Pipeline Development
* Data Processing and Statistical Analysis
* Scientific Data Visualization
* Linux-Based Workflow Management
* Reproducible Research Practices
* Primer Designing

## Repository Contents

* Workflow Documentation
* Pipeline Logic and Methodology
* Example Commands and Scripts
* Software Environment Information
* Computational Workflow Diagrams

## Repository Notice

This repository does not include raw sequencing data, processed datasets, experimental results, figures, or unpublished biological findings. It is intended to showcase the computational workflow, analytical methodology, and bioinformatics tools used during the study while maintaining research confidentiality.

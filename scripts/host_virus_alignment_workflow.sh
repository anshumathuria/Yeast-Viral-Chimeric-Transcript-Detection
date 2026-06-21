#!/bin/bash

# =============================================================================
# Automated RNA-seq Host–Virus Junction Detection Pipeline
# =============================================================================
#
# Project:
# Understanding the Crosstalk Between Viral Replication and Yeast Genome Ageing
#
# Description:
# This Bash workflow automates the complete RNA-seq processing pipeline for
# the detection of host–virus chimeric junctions. Raw paired-end sequencing
# reads undergo quality assessment using FastQC, alignment to the
# Saccharomyces cerevisiae reference genome using STAR, extraction of
# unmapped reads, remapping to viral genomes (L-A and M1), and host–virus
# chimeric alignment analysis using a combined reference genome.
#
# The pipeline generates STAR Chimeric.out.junction files containing
# candidate host–virus fusion events. These junction files are subsequently
# used for downstream analyses including junction filtering, PCR duplicate
# removal, fusion sequence reconstruction, primer design, and experimental
# validation.
#
# Workflow:
#   1. FastQC quality assessment
#   2. STAR alignment to yeast genome
#   3. Extraction of unmapped reads
#   4. Viral genome alignment (L-A and M1)
#   5. Host–virus chimeric alignment
#   6. Generation of Chimeric.out.junction files
#
# Input:
#   - Paired-end RNA-seq FASTQ files
#   - Yeast reference genome index
#   - Viral genome index
#   - Combined host–virus reference genome index
#
# Output:
#   - FastQC reports
#   - Yeast alignment BAM files
#   - Unmapped read files
#   - Viral alignment files
#   - STAR Chimeric.out.junction files
#
# Author: Anshu Mathuria
# MSc Biotechnology & Bioinformatics
# Institute of Bioinformatics and Applied Biotechnology (IBAB)
#
# =============================================================================
#!/bin/bash

# =====================================================
# RNA-seq Host-Virus Junction Detection Pipeline
# =====================================================

THREADS=16

mkdir -p fastqc
mkdir -p yeast_alignment
mkdir -p unmapped_reads
mkdir -p viral_alignment
mkdir -p chimeric_alignment

# =====================================================
# STEP 1: Run FastQC on all samples
# =====================================================

echo "Running FastQC..."

fastqc *.fastq.gz \
    -o fastqc \
    -t $THREADS

# =====================================================
# STEP 2–5: Process each sample
# =====================================================

for R1 in *_R1.fastq.gz
do

    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${SAMPLE}_R2.fastq.gz"

    echo "======================================="
    echo "Processing Sample: $SAMPLE"
    echo "======================================="

    # ----------------------------------
    # Step 2: Align to Yeast Genome
    # ----------------------------------

    STAR \
        --runThreadN $THREADS \
        --genomeDir yeast_index \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix yeast_alignment/${SAMPLE}_

    # ----------------------------------
    # Step 3: Extract Unmapped Reads
    # ----------------------------------

    samtools view \
        -b \
        -f 4 \
        yeast_alignment/${SAMPLE}_Aligned.sortedByCoord.out.bam \
        > unmapped_reads/${SAMPLE}_unmapped.bam

    samtools fastq \
        unmapped_reads/${SAMPLE}_unmapped.bam \
        > unmapped_reads/${SAMPLE}_unmapped.fastq

    # ----------------------------------
    # Step 4: Map Unmapped Reads to Virus
    # ----------------------------------

    STAR \
        --runThreadN $THREADS \
        --genomeDir viral_index \
        --readFilesIn unmapped_reads/${SAMPLE}_unmapped.fastq \
        --outFileNamePrefix viral_alignment/${SAMPLE}_

    # ----------------------------------
    # Step 5: Host-Virus Chimeric Mapping
    # ----------------------------------

    STAR \
        --runThreadN $THREADS \
        --genomeDir combined_host_virus_index \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions \
        --chimOutJunctionFormat 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix chimeric_alignment/${SAMPLE}_

    echo "Completed: $SAMPLE"

done

echo "======================================="
echo "All samples processed successfully"
echo "======================================="

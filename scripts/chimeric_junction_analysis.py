# =============================================================================
# Junction Analysis of Host–Virus Chimeric Reads
# =============================================================================
#
# Project:
# Understanding the Crosstalk Between Viral Replication and Yeast Genome Ageing
#
# Description:
# This script processes STAR chimeric junction output to identify high-confidence
# viral-host fusion events from RNA-seq data. The workflow filters viral-host
# junctions, removes potential PCR duplicates based on genomic breakpoints,
# selects confident candidate fusion sites, and reconstructs host-virus fusion
# sequences using BioPython and reference genomes. The generated fusion
# sequences can be used for downstream validation experiments such as RT-PCR,
# RT-qPCR primer design, and experimental confirmation of chimeric junctions.
#
# Biological System:
# Saccharomyces cerevisiae infected with L-A and M1 dsRNA viruses
#
# Experimental Conditions:
# Wild Type (WT) and gnc4 mutant strains under different nutrient conditions
# and time points.
#
# Author: Anshu Mathuria
# MSc Biotechnology & Bioinformatics
# Institute of Bioinformatics and Applied Biotechnology (IBAB)
#
# Purpose:
# Research project conducted as part of MSc thesis work to investigate
# potential interactions between viral replication and yeast genome ageing.
#
# =============================================================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
from Bio import SeqIO
from Bio.Seq import Seq

# ==============================================================================
# 1. Chimeric Junction Analysis & PCR Deduplication
# ==============================================================================

# Load the STAR chimeric junction output
cols = ['chr_donor', 'brkpt_donor', 'strand_donor', 'chr_acceptor', 'brkpt_acceptor', 'strand_acceptor', 'junction_type', 'read_name', 'first_base_pe', 'CIGAR']
df_junc = pd.read_csv('combined_chimeric_junctions.tsv', sep='\t', header=None, names=cols)

# Filter for Viral-Host junctions only (e.g., L-A virus to Yeast)
viral_host_mask = (df_junc['chr_donor'].str.contains('virus')) & (~df_junc['chr_acceptor'].str.contains('virus'))
df_viral_host = df_junc[viral_host_mask].copy()

# Remove PCR Duplicates based on exact breakpoints
df_dedup = df_viral_host.drop_duplicates(subset=['chr_donor', 'brkpt_donor', 'chr_acceptor', 'brkpt_acceptor'])

# Isolate the top 13 confident junctions
confident_junctions = df_dedup.head(13).copy()

# ==============================================================================
# 2. BioPython: Fusion Sequence Construction for Primer Design
# ==============================================================================

# Load reference genomes into memory as dictionaries using BioPython
print("Loading reference genomes into memory...")
yeast_genome = SeqIO.to_dict(SeqIO.parse("reference/yeast_genome.fasta", "fasta"))
viral_genome = SeqIO.to_dict(SeqIO.parse("reference/la_virus_genome.fasta", "fasta"))

def construct_fusion_sequence(row, padding=150):
    """
    Uses BioPython to fetch the sequence upstream of the viral donor breakpoint 
    and downstream of the yeast acceptor breakpoint, creating the exact in vivo fusion sequence.
    """
    donor_chr = row['chr_donor']
    acceptor_chr = row['chr_acceptor']
    donor_bp = int(row['brkpt_donor'])
    acceptor_bp = int(row['brkpt_acceptor'])

    try:
        # Extract viral sequence (upstream of junction)
        # Using max(0) prevents negative indexing if the breakpoint is near the chromosome start
        donor_seq = viral_genome[donor_chr].seq[max(0, donor_bp - padding) : donor_bp]
        
        # Extract yeast sequence (downstream of junction)
        acceptor_seq = yeast_genome[acceptor_chr].seq[acceptor_bp : acceptor_bp + padding]
        
        # Concatenate to form the fusion sequence
        fusion_seq = donor_seq + acceptor_seq
        
        # BioPython allows for easy reverse complementation if the junction is on the negative strand
        if row['strand_donor'] == '-':
             fusion_seq = fusion_seq.reverse_complement()
            
        return str(fusion_seq)
        
    except KeyError as e:
        return f"Error: Contig {e} not found in reference fasta."

print("Constructing fusion sequences for top 13 junctions...")
confident_junctions['fusion_seq'] = confident_junctions.apply(construct_fusion_sequence, axis=1)

# Save this to a CSV so you can take it to Primer3 or your preferred primer design tool
confident_junctions[['chr_donor', 'brkpt_donor', 'chr_acceptor', 'brkpt_acceptor', 'fusion_seq']].to_csv('rt_qpcr_primer_targets.csv', index=False)
print("Fusion sequences saved to 'rt_qpcr_primer_targets.csv' successfully.")

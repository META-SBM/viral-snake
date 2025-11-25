# Base Snakefile - Core pipeline configuration
# DO NOT EDIT - tracked in git

import os
from pathlib import Path

# ============================================================================
# Configuration
# ============================================================================

configfile: "config/config.yaml"

# Extract database paths
DATABASES = config['databases']
THREADS = config['threads']

# Default thread allocation
DEFAULT_THREADS = {
    'megahit': 32,
    'diamond': 32,
    'kraken2_contigs': 32,
    'kraken2': 16,
    'mmseqs2': 16,
    'blast': 6,
    'minimap2': 6,
    'cutadapt': 8,
    'quast': 8,
    'fastqc': 4,
    'seqkit': 4,
}

# Merge user config with defaults (user config overrides)
THREADS = {**DEFAULT_THREADS, **config.get('threads', {})}

# Set working directory
WORKDIR = config.get("workdir")
workdir: WORKDIR

# ============================================================================
# Include All Modules
# ============================================================================

# QC
include: "modules/qc/cutadapt/cutadapt.smk"
include: "modules/qc/cutadapt/cutadapt_barcode_rescue.smk"
include: "modules/qc/fastqc.smk"
include: "modules/qc/read_stats.smk"
include: "modules/qc/subsample.smk"

# Assembly
include: "modules/assembly/megahit.smk"
include: "modules/assembly/metaspades.smk"
include: "modules/assembly/qc.smk"
include: "modules/assembly/anvio.smk"

# Taxonomic classification
include: "modules/kraken2/kraken2.smk"
include: "modules/kraken2/kraken2_contigs.smk"

# Homology search
include: "modules/blast/blast.smk"
include: "modules/diamond/diamond.smk"
include: "modules/diamond/diamond_filter.smk"
include: "modules/diamond/diamond_lca_taxonkit.smk"

# Mapping and clustering
include: "modules/minimap2/minimap2.smk"
include: "modules/mmseqs2/cluster.smk"
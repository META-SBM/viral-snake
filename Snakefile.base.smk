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

# Set working directory
WORKDIR = config.get("workdir")
workdir: WORKDIR

# ============================================================================
# Include All Modules
# ============================================================================

# QC
include: "modules/qc/cutadapt.smk"
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
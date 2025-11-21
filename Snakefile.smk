# Snakefile for tick metavirome assembly pipeline

import os
import sys
from pathlib import Path

from viral_snake.utils import get_samples, discover_references
from viral_snake.collections import load_collections_from_dir, validate_collections

include: "Snakefile.base.smk"

# ============================================================================
# Target Rules
# ============================================================================

# Pipeline parameters
QC_FILTER = "raw__cutadapt_no_mgi_min_len_90"
QC_FILTER = "raw"
ASSEMBLERS = ["megahit"]
MIN_CONTIG_LENGTH = 800

# ============================================================================
# Discovery
# ============================================================================

# Discover samples
SAMPLES = get_samples(os.path.join(WORKDIR, 'reads/raw'))
print(f"Found {len(SAMPLES)} samples: {SAMPLES[:5]}..." if len(SAMPLES) > 5 else f"Found samples: {SAMPLES}")

# Discover collections
COLLECTION_DIR = Path(WORKDIR) / "sample_collections"
COLLECTIONS = load_collections_from_dir(COLLECTION_DIR)
print(f"Found {len(COLLECTIONS)} collections: {list(COLLECTIONS.keys())}")

if COLLECTIONS:
    validate_collections(COLLECTIONS, SAMPLES)

# Discover reference genomes
REFERENCE_GENOMES = discover_references(
    os.path.join(WORKDIR, 'refseq_reference'),
    pattern='*/ncbi_dataset/data/genomic.fna'
)
print(f"Found {len(REFERENCE_GENOMES)} reference genomes")


rule all:
    input:
       # # Subsample for test
       # expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
       #         qc_filter='raw__cutadapt_no_mgi_min_len_90__subsample_100000', sample=SAMPLES[0:10]),
        # QC reports
       expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
               qc_filter=QC_FILTER, sample=SAMPLES),
        
        # Co-assemblies with DIAMOND annotation
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/hits_with_taxonomy.tsv",
               assembler=ASSEMBLERS, 
               collection=COLLECTIONS, 
               min_len=MIN_CONTIG_LENGTH),
        
        # Filtered DIAMOND results
        expand("assembly/{assembler}/raw/{sample}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/{filter_preset}/hits.tsv",
               assembler=ASSEMBLERS,
               collection=COLLECTIONS,
               min_len=MIN_CONTIG_LENGTH,
               sample=SAMPLES,
               filter_preset=["viral_strict", "bacterial"]),
        
        # Reference mapping
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/alignments.sorted.bam",
               assembler=ASSEMBLERS,
               collection=COLLECTIONS,
               min_len=MIN_CONTIG_LENGTH,
               reference_org=REFERENCE_GENOMES),

        # LCA for all DIAMOND hits
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/contig_lca_detailed.tsv",
               assembler=ASSEMBLERS,
               collection=COLLECTIONS,
               min_len=MIN_CONTIG_LENGTH),
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
QC_FILTER = "raw__cutadapt_mgi_virome"
ASSEMBLERS = ["megahit"]
MIN_CONTIG_LENGTH = 700

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
               qc_filter='raw', sample=SAMPLES),

        expand('qc/read_stats/collections/ALL_RAW_read_stats.tsv')

                
        # expand(
        #     "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/download_summary.txt",
        #     prefix="co_assembly/megahit/ALL_SAMPLES_MERGED",
        #     min_len=700,
        #     preset="faster",
        #     database="NR",
        #     filter_preset="viral_strict"
        # ),

        # Reference mapping
        # expand(
        #     "alignment/strobealign__default/references/ICTV_DB/ICTV_DB.sanitized/__reads__/{query_qc}/{sample}/coverage.tsv",
        #     query_qc=QC_FILTER,
        #     sample=SAMPLES
        # ),

        # # Filtered DIAMOND results
        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/hits_with_taxonomy.tsv",
        #        assembler=ASSEMBLERS,
        #        qc_filter=QC_FILTER,
        #        min_len=MIN_CONTIG_LENGTH,
        #        sample=SAMPLES[110:120]),
        
    #    #  # Co-assemblies with DIAMOND annotation
    #     expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/hits_with_taxonomy.tsv",
    #            assembler=ASSEMBLERS, 
    #            collection=['ALL_SAMPLES_MERGED'], 
    #            min_len=MIN_CONTIG_LENGTH),
        
       #  # Filtered DIAMOND results
       #  expand("assembly/{assembler}/raw/{sample}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/{filter_preset}/hits.tsv",
       #         assembler=ASSEMBLERS,
       #         collection=COLLECTIONS,
       #         min_len=MIN_CONTIG_LENGTH,
       #         sample=SAMPLES,
       #         filter_preset=["viral_strict", "bacterial"]),
        
       #  # Reference mapping
       #  expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/alignments.sorted.bam",
       #         assembler=ASSEMBLERS,
       #         collection=COLLECTIONS,
       #         min_len=MIN_CONTIG_LENGTH,
       #         reference_org=REFERENCE_GENOMES),

       #  # LCA for all DIAMOND hits
       #  expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/contig_lca_detailed.tsv",
       #         assembler=ASSEMBLERS,
       #         collection=COLLECTIONS,
       #         min_len=MIN_CONTIG_LENGTH),
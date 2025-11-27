# Snakefile for tick metavirome assembly pipeline

import os
import sys
from pathlib import Path

from viral_snake.utils import get_samples, discover_references
from viral_snake.collections import load_collections_from_dir, validate_collections

include: "Snakefile.base.smk"



# Pipeline parameters
QC_FILTER = "raw__cutadapt_mgi_virome4"
QC_FILTER = "raw__barcode_rescue"
ASSEMBLERS = ["megahit"]
MIN_CONTIG_LENGTH = 700


# ============================================================================
# Adapters
# ============================================================================

# Base adapters
ADAPTER_R1_BASE = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2_BASE = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Load data
SAMPLE_SHEET = pd.read_csv("config/sample_sheet.tsv", sep='\t')
BARCODES = pd.read_csv("config/UDI_barcodes_reverse.csv", header=None)
BARCODES.columns = ['barcode_num', 'index']

# Create mappings
SAMPLE_TO_BARCODE_NUM = dict(zip(SAMPLE_SHEET['ID'], SAMPLE_SHEET['Sample_ID']))
BARCODE_NUM_TO_INDEX = dict(zip(BARCODES['barcode_num'], BARCODES['index']))

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


# ============================================================================
# Target Rules
# ============================================================================

rule all:
    input:
       # # Subsample for test
       # expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
       #         qc_filter='raw__cutadapt_no_mgi_min_len_90__subsample_100000', sample=SAMPLES[0:10]),
        # QC reports
        expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
               qc_filter=QC_FILTER, sample=SAMPLES),

        expand("reads/raw__barcode_rescue/{sample}_R1.fastq.gz", sample=SAMPLES),

        expand("qc/fastqc/{qc_filter}/{sample}_R{read}_fastqc.html",
               qc_filter='raw', sample=SAMPLES, read = ['1', '2']),
        expand("kraken2/{confidence}/{qc_filter}/{sample}.bracken",
                qc_filter=QC_FILTER, sample=SAMPLES, confidence = '0.5'
            ),
        # expand("feature_tables/{feature_table_id}/taxonomy_table.tsv",
        #     feature_table_id='bracken-species-all-0.5-min-len-90'),
        # 
        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/contig_summary.tsv",
        #        assembler=ASSEMBLERS,
        #        qc_filter=QC_FILTER,
        #        min_len=MIN_CONTIG_LENGTH,
        #        sample=SAMPLES),
        
        # Co-assemblies with DIAMOND annotation
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/hits_with_taxonomy.tsv",
               assembler=ASSEMBLERS, 
               collection=['ALL_SAMPLES_MERGED_RESCUED'], 
               min_len=MIN_CONTIG_LENGTH),

        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/viral_strict/hits.tsv",
               assembler=ASSEMBLERS, 
               collection=['ALL_SAMPLES_MERGED_RESCUED'], 
               min_len=MIN_CONTIG_LENGTH),
        
        expand(
            "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/download_summary.txt",
            prefix="co_assembly/megahit/ALL_SAMPLES_MERGED_RESCUED",
            min_len=700,
            preset="faster",
            database="NR",
            filter_preset="viral_strict"
        )
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
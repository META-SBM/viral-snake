# User Snakefile - Define your workflow targets here
# Copy this to "Snakefile" and customize for your project

from viral_snake.utils import get_samples, discover_references
from viral_snake.collections import load_collections_from_dir, validate_collections

# ============================================================================
# Include Base Pipeline
# ============================================================================

include: "Snakefile.base.smk"

# ============================================================================
# Sample and Reference Discovery
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
# Pipeline Parameters
# ============================================================================

QC_FILTER = "raw__cutadapt_no_mgi_min_len_90"
ASSEMBLERS = ["megahit"]
MIN_CONTIG_LENGTH = 800
DIAMOND_PRESET = "faster"

# ============================================================================
# Target Rule - Define what you want to generate
# ============================================================================

rule all:
    input:
        # Subsample for testing (optional - comment out for full run)
        expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
               qc_filter='raw__cutadapt_no_mgi_min_len_90__subsample_100000',
               sample=SAMPLES[0:10]),

        # QC reports for all samples
        expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
               qc_filter=QC_FILTER,
               sample=SAMPLES),
        
        # Co-assemblies with DIAMOND annotation
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_{preset}/NR/hits_with_taxonomy.tsv",
               assembler=ASSEMBLERS,
               collection=COLLECTIONS,
               min_len=MIN_CONTIG_LENGTH,
               preset=DIAMOND_PRESET),
        
        # Filtered DIAMOND results (viral and bacterial)
        expand("assembly/{assembler}/raw/{sample}/contigs_formatted_minlen_{min_len}/diamond_{preset}/NR/{filter_preset}/hits.tsv",
               assembler=ASSEMBLERS,
               sample=SAMPLES,
               min_len=MIN_CONTIG_LENGTH,
               preset=DIAMOND_PRESET,
               filter_preset=["viral_strict", "bacterial"]),
        
        # Reference genome mapping
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/alignments.sorted.bam",
               assembler=ASSEMBLERS,
               collection=COLLECTIONS,
               min_len=MIN_CONTIG_LENGTH,
               reference_org=REFERENCE_GENOMES),

        # LCA for DIAMOND hits
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_{preset}/NR/contig_lca_detailed.tsv",
               assembler=ASSEMBLERS,
               collection=COLLECTIONS,
               min_len=MIN_CONTIG_LENGTH,
               preset=DIAMOND_PRESET),
# Snakefile for tick metavirome assembly pipeline

import sys
from pathlib import Path

from viral_snake.utils import get_samples, discover_references
from viral_snake.collections import load_collections_from_dir, validate_collections

# Get workdir from config or use default
WORKDIR = config.get("workdir")
workdir: WORKDIR

# Discover samples
SAMPLES = get_samples(os.path.join(WORKDIR, 'reads/raw'))

print(f"Found {len(SAMPLES)} samples: {SAMPLES[:5]}..." if len(SAMPLES) > 5 else f"Found samples: {SAMPLES}")

# Discover collections
COLLECTION_DIR = Path(WORKDIR) / "sample_collections"
COLLECTIONS = load_collections_from_dir(COLLECTION_DIR)
print(f"Found {len(COLLECTIONS)} collections: {list(COLLECTIONS.keys())}")

# Validate collections against available samples
if COLLECTIONS:
    validate_collections(COLLECTIONS, SAMPLES)

# Discover reference genomes
refs = discover_references(
    os.path.join(WORKDIR, 'refseq_reference'),
    pattern='*/ncbi_dataset/data/genomic.fna'
)
print(f"Found {len(refs)} reference genomes")

# Pipeline configuration
QC_FILTER = "raw__cutadapt_no_mgi_min_len_90"
ASSEMBLERS = ["megahit"]

# Manual collection list (if needed)
COOLS = [
    'gus_khrustalny_district__novoopokino_village_dermacentor_reticulatus',
    'gvardeysky_district__konstantinovka_settlement_dermacentor_reticulatus',
    'chuguyevsky_district__zhuravlevka_river_ixodes_persulcatus',
    'chuguyevsky_district__zhuravlevka_river_haemophysalis_japonica',
    'blagoveshchensky_district__jsc_polief_dermacentor'
]

# Include all submodules
include: "modules/qc/cutadapt.smk"
include: "modules/qc/fastqc.smk"
include: "modules/qc/read_stats.smk"

include: "modules/assembly/megahit.smk"
include: "modules/assembly/metaspades.smk"
include: "modules/assembly/qc.smk"

include: "modules/kraken2/kraken2.smk"
include: "modules/kraken2/kraken2_contigs.smk"

include: "modules/blast/blast.smk"
include: "modules/diamond/diamond.smk"
include: "modules/minimap2.smk"

include: "modules/mmseqs2/cluster.smk"
include: "modules/assembly/anvio.smk"


rule all:
    input:
        # FastQC reports
        expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
               qc_filter=QC_FILTER, sample=SAMPLES),
        
        # Co-assemblies
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/hits_with_taxonomy.tsv",
               assembler=ASSEMBLERS, collection=COLLECTIONS, min_len=800),
        
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/alignments.sorted.bam",
               assembler=ASSEMBLERS, collection=COLLECTIONS, min_len=800, reference_org=refs)
import os
from pathlib import Path
from viral_snake.utils import get_samples
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()

include: "Snakefile.base.smk"

# ============================================================================
# Multi-Dataset Configuration
# ============================================================================

DATASETS = {
    # 'TEST': {
    #     'fs_prefix': '/mnt/mgx/DATASETS/INTERNAL/VIROME',
    #     'qc_filter': 'raw__cutadapt_mgi_virome4',
    #     'min_contig_length': 700,
    #     'assemblers': ['megahit'],
    #     'samples': []
    # },
    'RUN3': {
        'fs_prefix': '/mnt/mgx/DATASETS/INTERNAL/VIROME',
        'qc_filter': 'raw__cutadapt_mgi_virome',
        'min_contig_length': 700,
        'assemblers': ['megahit'],
        'samples': []  # Populated by discovery
    },
    'VIROME4': {
        'fs_prefix': '/mnt/mgx/DATASETS/INTERNAL/VIROME',
        'qc_filter': 'raw__cutadapt_mgi_virome4',
        'min_contig_length': 700,
        'assemblers': ['megahit'],
        'samples': []
    },
    'VIROME5': {
        'fs_prefix': '/mnt/mgx/DATASETS/INTERNAL/VIROME',
        'qc_filter': 'raw__cutadapt_mgi_virome4',
        'min_contig_length': 700,
        'assemblers': ['megahit'],
        'samples': []
    },
    'VIROME6': {
        'fs_prefix': '/mnt/mgx/DATASETS/INTERNAL/VIROME',
        'qc_filter': 'raw__cutadapt_mgi_virome',
        'min_contig_length': 700,
        'assemblers': ['megahit'],
        'samples': []  # Populated by discovery
    },
}

# ============================================================================
# Sample Discovery
# ============================================================================

console.print()
console.print(Panel.fit(
    "[bold cyan]Multi-Dataset Pipeline[/bold cyan]\n"
    f"[dim]Datasets:[/dim] {len(DATASETS)}",
    border_style="cyan"
))
console.print()

# Create table for dataset info
table = Table(show_header=True, header_style="bold cyan")
table.add_column("Dataset", style="cyan", width=12)
table.add_column("Samples", justify="right", style="green")
table.add_column("QC Filter", style="yellow")
table.add_column("Min Length", justify="right", style="magenta")

total_samples = 0

for dataset_name, config in DATASETS.items():
    fs_prefix = config['fs_prefix']
    sample_dir = f"{fs_prefix}/{dataset_name}/reads/raw"
    
    if not os.path.exists(sample_dir):
        console.print(f"[bold red]✗ Error:[/bold red] Sample directory not found: {sample_dir}")
        raise FileNotFoundError(f"Sample directory not found: {sample_dir}")
    
    samples = get_samples(sample_dir)
    config['samples'] = samples
    total_samples += len(samples)
    
    table.add_row(
        dataset_name,
        str(len(samples)),
        config['qc_filter'],
        str(config['min_contig_length'])
    )

console.print(table)
console.print(f"\n[bold]Total samples across all datasets:[/bold] [green]{total_samples}[/green]\n")

# ============================================================================
# Targets
# ============================================================================

# Build target list
target_list = []

# Raw read counts for all datasets
for dataset in DATASETS:
    fs_prefix = DATASETS[dataset]['fs_prefix']
    samples = DATASETS[dataset]['samples']
    
    for sample in samples:
        target_list.append(
            f"{fs_prefix}/{dataset}/qc/read_stats/raw/{sample}_read_counts.tsv"
        )


target_list = []

for dataset in DATASETS:
    fs_prefix = DATASETS[dataset]['fs_prefix']
    qc_filter = DATASETS[dataset]['qc_filter']
    samples = DATASETS[dataset]['samples']
    min_len = DATASETS[dataset]['min_contig_length']
    assembler = DATASETS[dataset]['assemblers'][0]
    
    for sample in samples:
        # Reference path: this sample's own assembly
        reference_path = f"assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}"
        
        # Kraken2 
        target_list.append(
            f"{fs_prefix}/{dataset}/kraken2/0.5/{qc_filter}/{sample}.bracken"
        )

        # # QC: Raw read counts
        # target_list.append(
        #     f"{fs_prefix}/{dataset}/qc/read_stats/raw/{sample}_read_counts.tsv"
        # )
        
        # # # Assembly: Individual contigs
        # target_list.append(
        #     f"{fs_prefix}/{dataset}/assembly/{assembler}/{qc_filter}/{sample}/contigs.fa"
        # )
        
        # Assembly: Formatted contigs
        target_list.append(
            f"{fs_prefix}/{dataset}/assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/contigs.fa"
        )
        
        # # Alignment: Map sample back to its own assembly
        # target_list.append(
        #     f"{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{qc_filter}/{sample}/alignments.sorted.bam"
        # )
        
        # Alignment: Coverage stats
        # target_list.append(
        #     f"{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{qc_filter}/{sample}/coverage.tsv"
        # )

console.print(f"[dim]Generated {len(target_list)} targets[/dim]\n")


minimap2_targets = []

reference_path = f"references/refseq/references.clean"

dataset = 'RUN3'
fs_prefix = DATASETS[dataset]['fs_prefix']
qc_filter = DATASETS[dataset]['qc_filter']
samples = DATASETS[dataset]['samples']
min_len = DATASETS[dataset]['min_contig_length']
assembler = DATASETS[dataset]['assemblers'][0]
for s in samples:
    query_path = f'assembly/megahit/{qc_filter}/{s}/contigs_formatted_minlen_700'
    # query_path = f'assembly/megahit/{qc_filter}/{s}'
    target = f"{fs_prefix}/{dataset}/alignment/minimap2_asm10/{reference_path}/__contigs__/{query_path}/alignments.paf"
    minimap2_targets.append(target)

    target = f"{fs_prefix}/{dataset}/alignment/minimap2_asm10/{reference_path}/__contigs__/{query_path}/reference_coverage.tsv"
    minimap2_targets.append(target)


reference_path = f"references/refseq/references.clean"
dataset = 'VIROME4'
fs_prefix = DATASETS[dataset]['fs_prefix']
qc_filter = DATASETS[dataset]['qc_filter']
samples = DATASETS[dataset]['samples']
min_len = DATASETS[dataset]['min_contig_length']
assembler = DATASETS[dataset]['assemblers'][0]
for s in samples:
    query_path = f'assembly/megahit/{qc_filter}/{s}/contigs_formatted_minlen_700'
    # query_path = f'assembly/megahit/{qc_filter}/{s}'
    target = f"{fs_prefix}/{dataset}/alignment/minimap2_asm10/{reference_path}/__contigs__/{query_path}/alignments.paf"
    minimap2_targets.append(target)

    target = f"{fs_prefix}/{dataset}/alignment/minimap2_asm10/{reference_path}/__contigs__/{query_path}/reference_coverage.tsv"
    minimap2_targets.append(target)

target_list.append(
    f"/mnt/mgx/DATASETS/INTERNAL/VIROME/RUN3/feature_tables/bracken-species-all/taxonomy_table.tsv"
)




res = expand('/mnt/mgx/DATASETS/INTERNAL/VIROME/{dataset}/qc/read_stats/collections/ALL_TRIMMED_DISCARD_read_stats.tsv', dataset=DATASETS.keys())
res = expand('/mnt/mgx/DATASETS/INTERNAL/VIROME/VIROME6/qc/read_stats/collections/{collections}_read_stats.tsv', collections=['ALL_RAW', 'ALL_TRIMMED_DISCARD'])
res = expand('/mnt/mgx/DATASETS/INTERNAL/VIROME/{dataset}/co_assembly/megahit/ALL_SAMPLES_MERGED/contigs_formatted_minlen_700/diamond_faster/NR/LCA.tsv', dataset=DATASETS.keys())

# res = expand("/mnt/mgx/DATASETS/INTERNAL/VIROME/VIROME6/feature_tables/bracken-species-all/heatmap_{heatmap_preset}.pdf", heatmap_preset=['viral_all', 'all_taxa_0.1'])

# Main rule
rule all:
    input:
        res[0]

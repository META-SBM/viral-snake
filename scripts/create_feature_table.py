#!/usr/bin/env python3
"""
Create feature table metadata file by discovering samples from reads directory.

Scans for fastq files and generates meta.yaml for abundance table creation.
"""

import click
import yaml
from pathlib import Path


@click.command()
@click.option(
    '--dataset-dir',
    '-d',
    required=True,
    help='Dataset base directory'
)
@click.option(
    '--qc-filter',
    '-q',
    required=True,
    help='QC filter name (e.g., raw__cutadapt)'
)
@click.option(
    '--feature-table-id',
    '-f',
    required=True,
    help='Feature table ID (e.g., species-all)'
)
@click.option(
    '--description',
    '-desc',
    required=True,
    help='Description of the feature table'
)
@click.option(
    '--abundance-metric',
    '-m',
    default='new_est_reads',
    help='Abundance metric to use (default: new_est_reads)'
)
@click.option(
    '--level',
    '-l',
    default='S',
    help='Taxonomic level (default: S for species)'
)
def create_feature_table(dataset_dir, qc_filter, feature_table_id, description, 
                        abundance_metric, level):
    """Create feature table metadata by discovering samples from reads."""
    
    dataset_path = Path(dataset_dir)
    
    # Find all R1 fastq files
    reads_dir = dataset_path / 'reads' / qc_filter
    r1_files = list(reads_dir.glob('*_R1.fastq.gz'))
    
    if not r1_files:
        raise FileNotFoundError(
            f"No *_R1.fastq.gz files found in {reads_dir}"
        )
    
    # Extract sample names: "KLP_168_R1.fastq.gz" -> "KLP_168"
    samples = sorted([f.name.replace('_R1.fastq.gz', '') for f in r1_files])
    
    # Create output directory inside dataset-dir
    output_dir = dataset_path / 'feature_tables' / f'bracken-{feature_table_id}'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create metadata
    meta = {
        'samples': samples,
        'qc_filter': qc_filter,
        'abundance_metric': abundance_metric,
        'level': level,
        'description': description
    }
    
    # Write meta.yaml
    meta_file = output_dir / 'meta.yaml'
    with open(meta_file, 'w') as f:
        yaml.dump(meta, f, default_flow_style=False, sort_keys=False)
    
    click.echo(f"✓ Created feature table metadata: {meta_file}")
    click.echo(f"  Samples: {len(samples)}")
    click.echo(f"  QC filter: {qc_filter}")
    click.echo(f"  Samples found: {', '.join(samples)}")


if __name__ == '__main__':
    create_feature_table()
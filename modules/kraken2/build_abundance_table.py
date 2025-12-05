#!/usr/bin/env python3
"""
Build abundance table from Bracken output files.

Takes multiple .bracken files and creates a single abundance table where:
- Rows: tax_id__{taxonomy_id}
- Columns: Sample names (extracted from filenames)
- Values: Specified abundance metric (default: new_est_reads)
"""

import click
import pandas as pd
from pathlib import Path


@click.command()
@click.option(
    '--input-files',
    '-i',
    multiple=True,
    required=True,
    help='Bracken output files (.bracken)'
)
@click.option(
    '--abundance-metric',
    '-m',
    default='new_est_reads',
    help='Column to extract from bracken files (default: new_est_reads)'
)
@click.option(
    '--output',
    '-o',
    required=True,
    help='Output TSV file path'
)
def build_abundance_table(input_files, abundance_metric, output):
    """Build abundance table from multiple Bracken files."""
    
    sample_data = {}
    
    # Process each bracken file
    for filepath in input_files:
        # Extract sample name: "path/to/sample.bracken" -> "sample"
        sample_name = Path(filepath).name.split('.')[0]
        
        # Read bracken file and extract taxonomy_id and metric column
        df = pd.read_csv(filepath, sep='\t')
        
        # Check if metric exists
        if abundance_metric not in df.columns:
            raise ValueError(
                f"Column '{abundance_metric}' not found in {filepath}. "
                f"Available columns: {', '.join(df.columns)}"
            )
        
        # Create Series with taxonomy_id as index
        series = df.set_index('taxonomy_id')[abundance_metric]
        sample_data[sample_name] = series
    
    # Combine all samples into single DataFrame
    # Missing taxa in some samples will be NaN
    abundance_df = pd.DataFrame(sample_data)
    
    # Fill missing values with 0
    abundance_df = abundance_df.fillna(0)
    
    # Convert to int if metric is read counts (not fractions)
    if abundance_metric in ['new_est_reads', 'added_reads', 'kraken_assigned_reads']:
        abundance_df = abundance_df.astype(int)
    
    # Format index with tax_id__ prefix
    abundance_df.index = 'tax_id__' + abundance_df.index.astype(str)
    
    # Write to TSV
    abundance_df.to_csv(output, sep='\t')
    
    click.echo(f"✓ Created abundance table: {output}")
    click.echo(f"  Samples: {len(sample_data)}")
    click.echo(f"  Taxa: {len(abundance_df)}")


if __name__ == '__main__':
    build_abundance_table()
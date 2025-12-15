#!/usr/bin/env python3
"""Plot mosdepth coverage results."""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import gzip
import pickle
import click
from pathlib import Path


def plot_mosdepth_coverage(
    regions_file,
    output_html,
    output_pickle,
    sample_name=None,
    verbose=True
):
    """
    Create coverage plots from mosdepth regions file.
    
    Args:
        regions_file: Path to mosdepth .regions.bed.gz file
        output_html: Path to save HTML plot
        output_pickle: Path to save pickle of figure
        sample_name: Sample name for title (extracted from path if None)
        verbose: Print progress messages
    
    Returns:
        dict: Statistics about coverage
    """
    if verbose:
        click.echo(f"Reading {regions_file}...", err=True)
    
    # Read the mosdepth regions file
    with gzip.open(regions_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t', names=['chrom', 'start', 'end', 'coverage'])
    
    # Get unique chromosomes/contigs
    contigs = df['chrom'].unique()
    n_contigs = len(contigs)
    
    if verbose:
        click.echo(f"Found {n_contigs} contig(s)", err=True)
        for contig in contigs:
            contig_df = df[df['chrom'] == contig]
            click.echo(f"  - {contig}: {len(contig_df):,} windows, "
                      f"length: {contig_df['end'].max():,} bp", err=True)
    
    # Extract sample name from path if not provided
    if sample_name is None:
        sample_name = Path(regions_file).parent.name
    
    # Create subplots - one row per contig
    fig = make_subplots(
        rows=n_contigs, cols=1,
        subplot_titles=[f"{contig}" for contig in contigs],
        vertical_spacing=0.08 / n_contigs if n_contigs > 1 else 0.1,
        shared_xaxes=False
    )
    
    # Statistics dictionary
    stats = {'contigs': {}, 'overall': {}}
    
    # Plot each contig separately
    for i, contig in enumerate(contigs, start=1):
        contig_df = df[df['chrom'] == contig].copy()
        mean_cov = contig_df['coverage'].mean()
        
        # Store stats
        stats['contigs'][contig] = {
            'length': int(contig_df['end'].max()),
            'windows': len(contig_df),
            'mean_coverage': float(mean_cov),
            'median_coverage': float(contig_df['coverage'].median()),
            'min_coverage': float(contig_df['coverage'].min()),
            'max_coverage': float(contig_df['coverage'].max()),
            'std_coverage': float(contig_df['coverage'].std())
        }
        
        # Add coverage trace
        fig.add_trace(
            go.Scatter(
                x=contig_df['start'],
                y=contig_df['coverage'],
                mode='lines',
                name=contig,
                line=dict(color='#3498db', width=1.5),
                fill='tozeroy',
                fillcolor='rgba(52, 152, 219, 0.2)',
                hovertemplate=f'<b>{contig}</b><br><b>Position:</b> %{{x:,}} bp<br><b>Coverage:</b> %{{y:.1f}}x<extra></extra>',
                showlegend=False
            ),
            row=i, col=1
        )
        
        # Add mean line for this contig
        fig.add_hline(
            y=mean_cov,
            line_dash="dash",
            line_color="rgba(231, 76, 60, 0.6)",
            line_width=1.5,
            annotation_text=f"Mean: {mean_cov:.1f}x",
            annotation_position="top right",
            annotation_font_size=10,
            row=i, col=1
        )
        
        # Update axes for this subplot
        fig.update_xaxes(title_text="Position (bp)", row=i, col=1, 
                        showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Coverage", row=i, col=1, 
                        showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    
    # Overall statistics
    stats['overall'] = {
        'total_windows': len(df),
        'mean_coverage': float(df['coverage'].mean()),
        'median_coverage': float(df['coverage'].median()),
        'min_coverage': float(df['coverage'].min()),
        'max_coverage': float(df['coverage'].max()),
        'std_coverage': float(df['coverage'].std())
    }
    
    # Update layout
    height = max(400, 300 * n_contigs)  # Dynamic height based on number of contigs
    fig.update_layout(
        height=height,
        width=1200,
        template='plotly_white',
        title_text=f"<b>{sample_name}</b> - Coverage by Contig",
        title_x=0.5,
        title_font_size=16,
        hovermode='x unified'
    )
    
    # Save HTML
    if verbose:
        click.echo(f"Saving HTML to {output_html}...", err=True)
    fig.write_html(output_html)
    
    # Save pickle
    if verbose:
        click.echo(f"Saving pickle to {output_pickle}...", err=True)
    with open(output_pickle, 'wb') as f:
        pickle.dump(fig, f)
    
    if verbose:
        click.echo("✓ Plot generation complete", err=True)
        click.echo(f"\nOverall statistics:", err=True)
        click.echo(f"  Mean coverage:   {stats['overall']['mean_coverage']:.2f}x", err=True)
        click.echo(f"  Median coverage: {stats['overall']['median_coverage']:.2f}x", err=True)
    
    return stats


@click.command()
@click.argument('regions_file', type=click.Path(exists=True))
@click.argument('output_html', type=click.Path())
@click.argument('output_pickle', type=click.Path())
@click.option('--sample-name', type=str, default=None,
              help='Sample name for plot title (default: extracted from path)')
@click.option('--quiet', is_flag=True,
              help='Suppress progress output')
def main(regions_file, output_html, output_pickle, sample_name, quiet):
    """
    Plot mosdepth coverage results.
    
    Takes a mosdepth regions.bed.gz file and creates interactive coverage plots
    with separate subplots for each contig/chromosome. Saves both HTML and pickle formats.
    
    Examples:
    
    \b
    # Basic usage
    plot_mosdepth.py coverage.regions.bed.gz plot.html plot.pickle
    
    \b
    # Custom sample name
    plot_mosdepth.py coverage.regions.bed.gz plot.html plot.pickle --sample-name "Sample_123"
    """
    plot_mosdepth_coverage(
        regions_file=regions_file,
        output_html=output_html,
        output_pickle=output_pickle,
        sample_name=sample_name,
        verbose=not quiet
    )


if __name__ == "__main__":
    main()
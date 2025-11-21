#!/usr/bin/env python3
"""Filter DIAMOND results by taxonomy and quality thresholds."""

import pandas as pd
import sys
import click


def filter_diamond_results(
    input_file,
    output_file,
    min_pident=50.0,
    min_length=100,
    max_evalue=1e-5,
    keep_domains=None,
    exclude_domains=None,
    min_bitscore=None,
    verbose=True
):
    """Core filtering logic."""
    df = pd.read_csv(input_file, sep="\t")
    
    if verbose:
        click.echo(f"Initial hits: {len(df)}", err=True)
    
    # Apply filters
    df = df[df['pident'] >= min_pident]
    if verbose:
        click.echo(f"After pident >= {min_pident}: {len(df)}", err=True)
    
    df = df[df['length'] >= min_length]
    if verbose:
        click.echo(f"After length >= {min_length}: {len(df)}", err=True)
    
    df = df[df['evalue'] <= max_evalue]
    if verbose:
        click.echo(f"After evalue <= {max_evalue}: {len(df)}", err=True)
    
    if min_bitscore is not None:
        df = df[df['bitscore'] >= min_bitscore]
        if verbose:
            click.echo(f"After bitscore >= {min_bitscore}: {len(df)}", err=True)
    
    # Taxonomic filters
    if keep_domains:
        df = df[df['Domain'].isin(keep_domains)]
        if verbose:
            click.echo(f"After keeping domains {keep_domains}: {len(df)}", err=True)
    
    if exclude_domains:
        df = df[~df['Domain'].isin(exclude_domains)]
        if verbose:
            click.echo(f"After excluding domains {exclude_domains}: {len(df)}", err=True)
    
    # Save
    df.to_csv(output_file, sep="\t", index=False)
    
    if verbose:
        click.echo(f"✓ Final hits: {len(df)}", err=True)
        click.echo(f"✓ Saved to: {output_file}", err=True)
    
    return len(df)


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
@click.option('--min-pident', default=50.0, type=float,
              help='Minimum percent identity (default: 50)')
@click.option('--min-length', default=100, type=int,
              help='Minimum alignment length (default: 100)')
@click.option('--max-evalue', default=1e-5, type=float,
              help='Maximum e-value (default: 1e-5)')
@click.option('--min-bitscore', type=float,
              help='Minimum bitscore')
@click.option('--keep-domains', multiple=True,
              help='Keep only these domains (can use multiple times)')
@click.option('--exclude-domains', multiple=True,
              help='Exclude these domains (can use multiple times)')
@click.option('--quiet', is_flag=True,
              help='Suppress progress output')
def main(input_file, output_file, min_pident, min_length, max_evalue,
         min_bitscore, keep_domains, exclude_domains, quiet):
    """
    Filter DIAMOND results by taxonomy and quality thresholds.
    
    Examples:
    
    \b
    # Keep only viral hits with high identity
    filter_diamond.py hits.tsv viral.tsv --keep-domains Viruses --min-pident 70
    
    \b
    # Exclude eukaryotic contamination
    filter_diamond.py hits.tsv filtered.tsv --exclude-domains Eukaryota
    
    \b
    # Multiple domains
    filter_diamond.py hits.tsv filtered.tsv \\
        --keep-domains Viruses \\
        --keep-domains Bacteria
    """
    keep_domains = list(keep_domains) if keep_domains else None
    exclude_domains = list(exclude_domains) if exclude_domains else None
    
    filter_diamond_results(
        input_file=input_file,
        output_file=output_file,
        min_pident=min_pident,
        min_length=min_length,
        max_evalue=max_evalue,
        keep_domains=keep_domains,
        exclude_domains=exclude_domains,
        min_bitscore=min_bitscore,
        verbose=not quiet
    )


if __name__ == "__main__":
    main()
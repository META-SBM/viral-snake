"""Feature table creation command"""

import click
from rich.panel import Panel
from rich.prompt import Confirm
from pathlib import Path
import yaml
from datetime import datetime

from . import (
    console,
    create_group,
    resolve_samples,
    apply_exclusions,
    validate_samples_exist,
    display_yaml_preview
)
from ... import __version__


def _expand_level(level):
    """Expand single-letter taxonomic level to full name"""
    levels = {
        'D': 'Domain',
        'P': 'Phylum',
        'C': 'Class',
        'O': 'Order',
        'F': 'Family',
        'G': 'Genus',
        'S': 'Species'
    }
    return levels.get(level.upper(), level)


@create_group.command(name='feature-table')
@click.argument('dataset_root', type=click.Path(exists=True))
@click.option('--table-id', '-t', required=True, 
              help='Feature table ID (e.g., species-all)')
@click.option('--qc-filter', '-q', required=True, 
              help='QC filter name (e.g., raw__cutadapt_mgi_virome)')
@click.option('--description', '-d', help='Feature table description (optional)')
@click.option('--level', '-l', default='S', 
              type=click.Choice(['D', 'P', 'C', 'O', 'F', 'G', 'S'], case_sensitive=False),
              help='Taxonomic level (default: S for species)')
@click.option('--abundance-metric', '-m', default='new_est_reads',
              type=click.Choice(['new_est_reads', 'fraction_total_reads'], case_sensitive=False),
              help='Abundance metric from Bracken (default: new_est_reads)')
@click.option('--samples', '-s', multiple=True, 
              help='Sample names (use ALL for all samples)')
@click.option('--sample-file', type=click.Path(exists=True), 
              help='File with sample names (one per line)')
@click.option('--exclude', '-x', multiple=True,
              help='Samples to exclude')
@click.option('--exclude-file', type=click.Path(exists=True),
              help='File with samples to exclude')
@click.option('--force', '-f', is_flag=True, help='Overwrite existing')
@click.option('--dry-run', is_flag=True, help='Preview without creating')
@click.option('--no-validate', is_flag=True, help='Skip validation')
def create_feature_table(dataset_root, table_id, qc_filter, description, level,
                        abundance_metric, samples, sample_file, exclude, exclude_file,
                        force, dry_run, no_validate):
    """Create feature table metadata for abundance analysis
    
    Discovers samples from reads directory and generates meta.yaml for 
    creating abundance matrices from Kraken2/Bracken results.
    
    \b
    Taxonomic levels:
      D = Domain (Bacteria, Viruses, etc.)
      P = Phylum
      C = Class
      O = Order
      F = Family
      G = Genus
      S = Species (default)
    
    \b
    Abundance metrics:
      new_est_reads        = Bracken-estimated read counts (default)
      fraction_total_reads = Fraction of total reads
    
    \b
    Examples:
    
    \b
    # Species-level table for all samples, excluding controls
    viral-snake create feature-table /path/to/dataset \\
        --table-id species-all \\
        --qc-filter raw__cutadapt_mgi_virome \\
        --description "Species-level abundance" \\
        --samples ALL \\
        --exclude negative_control --exclude positive_control
    
    \b
    # Genus-level from file with custom metric
    viral-snake create feature-table /path/to/dataset \\
        --table-id genus-fraction \\
        --qc-filter raw \\
        --level G \\
        --abundance-metric fraction_total_reads \\
        --sample-file samples.txt
    """
    
    dataset_path = Path(dataset_root)
    
    # Print header
    console.print()
    console.print(Panel.fit(
        f"[bold cyan]Creating Feature Table Metadata[/bold cyan]\n"
        f"[dim]Dataset:[/dim] {dataset_root}\n"
        f"[dim]Table ID:[/dim] {table_id}\n"
        f"[dim]Level:[/dim] {level} ({_expand_level(level)})\n"
        f"[dim]Metric:[/dim] {abundance_metric}",
        border_style="cyan"
    ))
    console.print()
    
    # ========================================================================
    # Resolve and validate samples
    # ========================================================================
    console.print("[bold]Sample Discovery[/bold]")
    console.print(f"  QC Filter: [cyan]{qc_filter}[/cyan]")
    console.print()
    
    try:
        # Resolve samples
        sample_list = resolve_samples(dataset_path, qc_filter, samples, sample_file)
        initial_count = len(sample_list)
        
        # Apply exclusions
        if exclude or exclude_file:
            sample_list, excluded_samples, not_found = apply_exclusions(
                sample_list, exclude, exclude_file
            )
            
            # Check we still have samples
            if not sample_list:
                raise click.ClickException(
                    "No samples remaining after exclusions!\n"
                    f"Started with {initial_count}, excluded {len(excluded_samples)}"
                )
        else:
            excluded_samples = []
        
        # Validate samples exist (unless skipped)
        if not no_validate:
            reads_dir = dataset_path / 'reads' / qc_filter
            valid_samples, missing_samples = validate_samples_exist(
                reads_dir, sample_list, show_progress=True
            )
            
            if missing_samples:
                console.print(f"[yellow]⚠ Warning:[/yellow] {len(missing_samples)} samples not found:")
                for s in missing_samples[:5]:
                    console.print(f"    [dim]• {s}[/dim]")
                if len(missing_samples) > 5:
                    console.print(f"    [dim]... and {len(missing_samples) - 5} more[/dim]")
                
                if not Confirm.ask(f"\n[yellow]Continue with {len(valid_samples)} valid samples?[/yellow]"):
                    raise click.Abort()
                
                sample_list = valid_samples
            else:
                console.print(f"[green]✓[/green] All {len(sample_list)} samples validated")
        
        # Show sample preview
        console.print(f"\n[dim]Final sample list (showing first 10):[/dim]")
        for i, s in enumerate(sorted(sample_list)[:10], 1):
            console.print(f"  [dim]{i:3d}.[/dim] {s}")
        if len(sample_list) > 10:
            console.print(f"  [dim]... and {len(sample_list) - 10} more[/dim]")
        console.print()
        
    except click.ClickException as e:
        console.print(f"[bold red]✗ Error:[/bold red] {e}")
        raise click.Abort()
    
    # ========================================================================
    # Create metadata structure
    # ========================================================================
    meta = {
        'samples': sorted(sample_list),
        'qc_filter': qc_filter,
        'abundance_metric': abundance_metric,
        'level': level
    }
    
    # Add optional fields
    if description:
        meta['description'] = description
    
    meta['created'] = datetime.now().isoformat()
    meta['created_by'] = f'viral-snake v{__version__}'
    
    # Add exclusion info if used
    if excluded_samples:
        meta['excluded_samples'] = sorted(excluded_samples)
        meta['samples_before_exclusion'] = initial_count
    
    # ========================================================================
    # Show preview
    # ========================================================================
    console.print("[bold]Preview:[/bold]\n")
    display_yaml_preview(meta, title=f"Feature Table: {table_id}")
    console.print()
    
    # Determine output path
    output_dir = dataset_path / 'feature_tables' / f'bracken-{table_id}'
    output_dir.mkdir(parents=True, exist_ok=True)
    meta_file = output_dir / 'meta.yaml'
    
    # ========================================================================
    # Dry run exit
    # ========================================================================
    if dry_run:
        console.print(f"[yellow]🔍 Dry run - would create:[/yellow] {meta_file}")
        console.print("[dim]Remove --dry-run to actually create the file[/dim]\n")
        return
    
    # ========================================================================
    # Check for existing file
    # ========================================================================
    if meta_file.exists() and not force:
        console.print(f"[yellow]⚠ File already exists:[/yellow] {meta_file}")
        if not Confirm.ask("[yellow]Overwrite?[/yellow]"):
            console.print("[dim]Aborted[/dim]\n")
            raise click.Abort()
    
    # ========================================================================
    # Write file
    # ========================================================================
    with open(meta_file, 'w') as f:
        yaml.dump(meta, f, default_flow_style=False, sort_keys=False)
    
    # ========================================================================
    # Success summary
    # ========================================================================
    summary_text = (
        f"[bold green]✓ Feature Table Metadata Created[/bold green]\n\n"
        f"[cyan]File:[/cyan] {meta_file}\n"
        f"[cyan]Table ID:[/cyan] {table_id}\n"
        f"[cyan]Level:[/cyan] {level} ({_expand_level(level)})\n"
        f"[cyan]Metric:[/cyan] {abundance_metric}\n"
        f"[cyan]Samples:[/cyan] {len(sample_list)}\n"
        f"[cyan]QC Filter:[/cyan] {qc_filter}"
    )
    
    if excluded_samples:
        summary_text += f"\n[yellow]Excluded:[/yellow] {len(excluded_samples)} samples"
    
    console.print(Panel.fit(summary_text, border_style="green"))
    
    # ========================================================================
    # Next steps
    # ========================================================================
    console.print("\n[bold]Next steps:[/bold]")
    console.print(f"  1. Run Kraken2: [cyan]snakemake kraken2/{qc_filter}/{{sample}}_report.txt[/cyan]")
    console.print(f"  2. Run Bracken: [cyan]snakemake bracken/{qc_filter}/{{sample}}_bracken_{level}.txt[/cyan]")
    console.print(f"  3. Create table: [cyan]snakemake feature_tables/bracken-{table_id}/abundance_table.tsv[/cyan]")
    console.print()
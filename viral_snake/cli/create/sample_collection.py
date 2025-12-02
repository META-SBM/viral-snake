"""Sample collection creation command"""

import click
from rich.panel import Panel
from rich.prompt import Confirm
from pathlib import Path
import yaml
from datetime import datetime

from . import (
    console, 
    create_group,
    validate_collection_name,
    resolve_samples,
    apply_exclusions,
    validate_samples_exist,
    display_yaml_preview
)
from ... import __version__


@create_group.command(name='sample-collection')
@click.argument('dataset_root', type=click.Path(exists=True))
@click.option('--name', '-n', required=True, help='Collection name')
@click.option('--qc-filter', '-q', required=True, help='QC filter to use')
@click.option('--description', '-d', help='Collection description (optional)')
@click.option('--samples', '-s', multiple=True, help='Sample names (use ALL for all)')
@click.option('--sample-file', type=click.Path(exists=True), help='File with sample names')
@click.option('--exclude', '-x', multiple=True, help='Samples to exclude')
@click.option('--exclude-file', type=click.Path(exists=True), help='File with exclusions')
@click.option('--force', '-f', is_flag=True, help='Overwrite existing')
@click.option('--dry-run', is_flag=True, help='Preview without creating')
@click.option('--no-validate', is_flag=True, help='Skip validation')
def create_sample_collection(dataset_root, name, qc_filter, description, samples,
                             sample_file, exclude, exclude_file, force, dry_run,
                             no_validate):
    """Create a sample collection for co-assembly
    
    Sample collections define groups of samples that will be co-assembled together.
    The collection can contain 1 or more samples.
    
    \b
    Examples:
    
    \b
    # All samples from a QC filter
    viral-snake create sample-collection /path/to/dataset \\
        --name tick_pools \\
        --qc-filter raw__cutadapt_mgi_virome \\
        --description "Tick samples for co-assembly" \\
        --samples ALL
    
    \b
    # Specific samples with exclusions
    viral-snake create sample-collection /path/to/dataset \\
        --name forest_sites \\
        --qc-filter raw \\
        --samples sample1 --samples sample2 --samples sample3 \\
        --exclude control1 --exclude control2
    
    \b
    # From file
    viral-snake create sample-collection /path/to/dataset \\
        --name selected_samples \\
        --qc-filter raw \\
        --sample-file samples.txt
    """
    
    dataset_path = Path(dataset_root)
    
    # Validate collection name
    try:
        validate_collection_name(name)
    except click.ClickException as e:
        console.print(f"[bold red]✗ Error:[/bold red] {e}")
        raise click.Abort()
    
    # Print header
    console.print()
    console.print(Panel.fit(
        f"[bold cyan]Creating Sample Collection[/bold cyan]\n"
        f"[dim]Dataset:[/dim] {dataset_root}\n"
        f"[dim]Name:[/dim] {name}\n"
        f"[dim]QC Filter:[/dim] {qc_filter}",
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
            not_found = []
        
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
        
        # Final check
        if not sample_list:
            raise click.ClickException("No valid samples remaining!")
        
        # Show sample preview
        console.print(f"\n[dim]Sample list (showing first 10):[/dim]")
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
    collection = {
        'name': name,
        'qc_filter': qc_filter,
        'samples': sorted(sample_list)
    }
    
    # Add optional fields
    if description:
        collection['description'] = description
    
    collection['created'] = datetime.now().isoformat()
    collection['created_by'] = f'viral-snake v{__version__}'
    
    # Add exclusion info if used
    if excluded_samples:
        collection['excluded_samples'] = sorted(excluded_samples)
        collection['samples_before_exclusion'] = initial_count
    
    # ========================================================================
    # Show preview
    # ========================================================================
    console.print("[bold]Preview:[/bold]\n")
    display_yaml_preview(collection, title=f"Sample Collection: {name}")
    console.print()
    
    # Determine output path
    config_dir = dataset_path / 'config' / 'sample_collections'
    config_dir.mkdir(parents=True, exist_ok=True)
    config_file = config_dir / f'{name}.yaml'
    
    # ========================================================================
    # Dry run exit
    # ========================================================================
    if dry_run:
        console.print(f"[yellow]🔍 Dry run - would create:[/yellow] {config_file}")
        console.print("[dim]Remove --dry-run to actually create the file[/dim]\n")
        return
    
    # ========================================================================
    # Check for existing file
    # ========================================================================
    if config_file.exists() and not force:
        console.print(f"[yellow]⚠ File already exists:[/yellow] {config_file}")
        if not Confirm.ask("[yellow]Overwrite?[/yellow]"):
            console.print("[dim]Aborted[/dim]\n")
            raise click.Abort()
    
    # ========================================================================
    # Write file
    # ========================================================================
    with open(config_file, 'w') as f:
        yaml.dump(collection, f, default_flow_style=False, sort_keys=False)
    
    # ========================================================================
    # Success summary
    # ========================================================================
    summary_text = (
        f"[bold green]✓ Sample Collection Created[/bold green]\n\n"
        f"[cyan]File:[/cyan] {config_file}\n"
        f"[cyan]Name:[/cyan] {name}\n"
        f"[cyan]QC Filter:[/cyan] {qc_filter}\n"
        f"[cyan]Samples:[/cyan] {len(sample_list)}"
    )
    
    if excluded_samples:
        summary_text += f"\n[yellow]Excluded:[/yellow] {len(excluded_samples)} samples"
    
    console.print(Panel.fit(summary_text, border_style="green"))
    
    # ========================================================================
    # Next steps
    # ========================================================================
    console.print("\n[bold]Next steps:[/bold]")
    console.print(f"  1. Run co-assembly: [cyan]snakemake co_assembly/megahit/{name}/contigs.fa[/cyan]")
    console.print(f"  2. Use in assembly collection: [cyan]--co-assembly-collections {name}[/cyan]")
    console.print()
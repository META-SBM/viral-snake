"""Assembly collection creation command"""

import click
from rich.panel import Panel
from rich.prompt import Confirm
from rich.progress import Progress, SpinnerColumn, TextColumn
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


@create_group.command(name='assembly-collection')
@click.argument('dataset_root', type=click.Path(exists=True))
@click.option('--name', '-n', required=True, help='Assembly collection name')
@click.option('--description', '-d', help='Collection description (optional)')
# Individual assemblies
@click.option('--assembler', default='megahit', 
              type=click.Choice(['megahit', 'metaspades']),
              help='Assembler to use (default: megahit)')
@click.option('--qc-filter', '-q', help='QC filter for individual assemblies')
@click.option('--min-len', type=int, default=800, 
              help='Minimum contig length (default: 800)')
@click.option('--samples', '-s', multiple=True, 
              help='Sample names (use ALL for all samples)')
@click.option('--sample-file', type=click.Path(exists=True), 
              help='File with sample names (one per line)')
@click.option('--exclude', '-x', multiple=True,
              help='Samples to exclude')
@click.option('--exclude-file', type=click.Path(exists=True),
              help='File with samples to exclude')
# Co-assemblies
@click.option('--co-assembly-collections', '-c', multiple=True, 
              help='Sample collection names for co-assemblies')
# Behavior
@click.option('--force', '-f', is_flag=True, help='Overwrite existing')
@click.option('--dry-run', is_flag=True, help='Preview without creating')
@click.option('--no-validate', is_flag=True, help='Skip validation')
def create_assembly_collection(dataset_root, name, description, assembler, 
                               qc_filter, min_len, samples, sample_file,
                               exclude, exclude_file,
                               co_assembly_collections, force, dry_run, 
                               no_validate):
    """Create an assembly collection configuration file
    
    The metadata file contains FULLY RESOLVED sample lists - no wildcards.
    Use ALL (uppercase) to auto-discover all samples.
    
    \b
    Examples:
    
    \b
    # All individual assemblies from raw QC filter
    viral-snake create assembly-collection /path/to/dataset \\
        --name all_ticks_raw \\
        --qc-filter raw \\
        --samples ALL
    
    \b
    # Specific samples with co-assemblies
    viral-snake create assembly-collection /path/to/dataset \\
        --name complete_dataset \\
        --description "All tick data" \\
        --qc-filter raw \\
        --samples ALL \\
        --co-assembly-collections tick_pools forest_pools
    
    \b
    # Just co-assemblies (no individual)
    viral-snake create assembly-collection /path/to/dataset \\
        --name coassemblies_only \\
        --co-assembly-collections tick_pools forest_pools
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
        f"[bold cyan]Creating Assembly Collection[/bold cyan]\n"
        f"[dim]Dataset:[/dim] {dataset_root}\n"
        f"[dim]Name:[/dim] {name}",
        border_style="cyan"
    ))
    console.print()
    
    sources = []
    stats = {'total_individual_assemblies': 0, 'total_co_assemblies': 0}
    
    # ========================================================================
    # Handle individual assemblies
    # ========================================================================
    if qc_filter:
        console.print("[bold]Individual Assemblies[/bold]")
        console.print(f"  Assembler: [cyan]{assembler}[/cyan]")
        console.print(f"  QC Filter: [cyan]{qc_filter}[/cyan]")
        console.print(f"  Min Length: [cyan]{min_len}[/cyan]")
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
                valid_samples, missing_samples = validate_samples_exist(reads_dir, sample_list)
                
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
            
            sources.append({
                'type': 'individual',
                'assembler': assembler,
                'qc_filter': qc_filter,
                'min_len': min_len,
                'samples': sorted(sample_list)
            })
            stats['total_individual_assemblies'] = len(sample_list)
            
            # Show sample preview
            console.print(f"\n[dim]Sample preview (showing first 10):[/dim]")
            for i, s in enumerate(sample_list[:10], 1):
                console.print(f"  [dim]{i:3d}.[/dim] {s}")
            if len(sample_list) > 10:
                console.print(f"  [dim]... and {len(sample_list) - 10} more[/dim]")
            console.print()
            
        except click.ClickException as e:
            console.print(f"[bold red]✗ Error:[/bold red] {e}")
            raise click.Abort()
    
    # ========================================================================
    # Handle co-assemblies
    # ========================================================================
    if co_assembly_collections:
        console.print("[bold]Co-assemblies[/bold]")
        console.print(f"  Assembler: [cyan]{assembler}[/cyan]")
        console.print(f"  Collections: [cyan]{', '.join(co_assembly_collections)}[/cyan]")
        console.print(f"  Min Length: [cyan]{min_len}[/cyan]")
        console.print()
        
        # Validate collections exist (unless skipped)
        if not no_validate:
            config_dir = dataset_path / 'config' / 'sample_collections'
            missing_collections = []
            for coll in co_assembly_collections:
                coll_file = config_dir / f"{coll}.yaml"
                if not coll_file.exists():
                    missing_collections.append(coll)
            
            if missing_collections:
                console.print(f"[bold red]✗ Error:[/bold red] Collections not found:")
                for c in missing_collections:
                    console.print(f"    • {c}")
                console.print(f"\n[dim]Create them first with:[/dim]")
                console.print(f"[dim]  viral-snake create sample-collection ...[/dim]")
                raise click.Abort()
        
        sources.append({
            'type': 'co_assembly',
            'assembler': assembler,
            'collections': sorted(list(co_assembly_collections)),
            'min_len': min_len
        })
        stats['total_co_assemblies'] = len(co_assembly_collections)
        console.print(f"[green]✓[/green] Added {len(co_assembly_collections)} co-assemblies\n")
    
    # ========================================================================
    # Validate we have sources
    # ========================================================================
    if not sources:
        console.print(
            "[bold red]✗ Error:[/bold red] No sources specified\n\n"
            "[yellow]You must specify at least one of:[/yellow]\n"
            "  • Individual assemblies: [cyan]--qc-filter[/cyan] and [cyan]--samples[/cyan]\n"
            "  • Co-assemblies: [cyan]--co-assembly-collections[/cyan]"
        )
        raise click.Abort()
    
    # ========================================================================
    # Create metadata structure
    # ========================================================================
    collection = {
        'name': name,
        'sources': sources,
        'stats': stats
    }
    
    # Add optional fields
    if description:
        collection['description'] = description
    
    collection['created'] = datetime.now().isoformat()
    collection['created_by'] = f'viral-snake v{__version__}'
    
    # ========================================================================
    # Show preview
    # ========================================================================
    console.print("[bold]Preview:[/bold]\n")
    display_yaml_preview(collection, title=f"Assembly Collection: {name}")
    console.print()
    
    # Determine output path
    config_dir = dataset_path / 'config' / 'assembly_collections'
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
        f"[bold green]✓ Assembly Collection Created[/bold green]\n\n"
        f"[cyan]File:[/cyan] {config_file}\n"
        f"[cyan]Name:[/cyan] {name}\n"
        f"[cyan]Individual assemblies:[/cyan] {stats['total_individual_assemblies']}\n"
        f"[cyan]Co-assemblies:[/cyan] {stats['total_co_assemblies']}"
    )
    
    console.print(Panel.fit(summary_text, border_style="green"))
    
    # ========================================================================
    # Next steps
    # ========================================================================
    console.print("\n[bold]Next steps:[/bold]")
    console.print(f"  1. Prepare collection: [cyan]snakemake assembly_collection/{name}/contigs.fa[/cyan]")
    console.print(f"  2. Run DIAMOND batch: [cyan]snakemake assembly_collection/{name}/diamond_faster/NR/hits.txt[/cyan]")
    console.print()
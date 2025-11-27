"""Commands for creating metadata and configuration files"""

import click
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.prompt import Confirm
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.syntax import Syntax
from pathlib import Path
import yaml
from datetime import datetime

from ..utils import get_samples
from .. import __version__

console = Console()


@click.group(name='create')
def create_group():
    """📝 Create metadata and configuration files"""
    pass


def display_yaml_preview(data, title="Preview"):
    """Display YAML data with syntax highlighting"""
    yaml_str = yaml.dump(data, default_flow_style=False, sort_keys=False)
    syntax = Syntax(yaml_str, "yaml", theme="monokai", line_numbers=True)
    console.print(Panel(syntax, title=f"📄 {title}", border_style="cyan"))


def resolve_samples(dataset_path, qc_filter, samples, sample_file):
    """Resolve sample list from various input methods
    
    Returns:
        list: Resolved sample names
        
    Raises:
        click.ClickException: If samples can't be resolved
    """
    reads_dir = dataset_path / 'reads' / qc_filter
    
    if not reads_dir.exists():
        raise click.ClickException(
            f"QC filter directory not found: {reads_dir}\n"
            f"Available filters: {', '.join([d.name for d in (dataset_path / 'reads').iterdir() if d.is_dir()])}"
        )
    
    # Handle different input methods
    if 'ALL' in samples or 'all' in samples:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console
        ) as progress:
            progress.add_task(description="Discovering all samples...", total=None)
            try:
                sample_list = get_samples(reads_dir)
            except FileNotFoundError as e:
                raise click.ClickException(str(e))
        
        console.print(f"[green]✓[/green] Discovered {len(sample_list)} samples")
        
    elif sample_file:
        console.print(f"[cyan]Reading samples from:[/cyan] {sample_file}")
        with open(sample_file) as f:
            sample_list = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        console.print(f"[green]✓[/green] Loaded {len(sample_list)} samples from file")
        
    else:
        sample_list = list(samples)
        console.print(f"[cyan]Using {len(sample_list)} specified samples[/cyan]")
    
    if not sample_list:
        raise click.ClickException(
            "No samples specified. Use:\n"
            "  --samples ALL (discover all)\n"
            "  --samples sample1 --samples sample2 (specific)\n"
            "  --sample-file samples.txt (from file)"
        )
    
    return sample_list


def validate_samples_exist(reads_dir, sample_list, show_progress=True):
    """Validate that all samples have R1 files
    
    Returns:
        tuple: (valid_samples, missing_samples)
    """
    valid = []
    missing = []
    
    if show_progress:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console
        ) as progress:
            task = progress.add_task(
                description=f"Validating {len(sample_list)} samples...", 
                total=len(sample_list)
            )
            
            for sample in sample_list:
                r1 = reads_dir / f"{sample}_R1.fastq.gz"
                if r1.exists():
                    valid.append(sample)
                else:
                    missing.append(sample)
                progress.advance(task)
    else:
        for sample in sample_list:
            r1 = reads_dir / f"{sample}_R1.fastq.gz"
            if r1.exists():
                valid.append(sample)
            else:
                missing.append(sample)
    
    return valid, missing


@create_group.command(name='assembly-collection')
@click.argument('dataset_root', type=click.Path(exists=True))
@click.option('--name', '-n', required=True, help='Assembly collection name')
@click.option('--description', '-d', help='Collection description')
# Individual assemblies
@click.option('--assembler', default='megahit', type=click.Choice(['megahit', 'metaspades']),
              help='Assembler to use')
@click.option('--qc-filter', '-q', help='QC filter for individual assemblies')
@click.option('--min-len', type=int, default=800, help='Minimum contig length')
@click.option('--samples', '-s', multiple=True, 
              help='Sample names (use ALL for all samples)')
@click.option('--sample-file', type=click.Path(exists=True), 
              help='File with sample names (one per line)')
# Co-assemblies
@click.option('--co-assembly-collections', '-c', multiple=True, 
              help='Sample collection names for co-assemblies')
# Behavior
@click.option('--force', '-f', is_flag=True, help='Overwrite existing file without asking')
@click.option('--dry-run', is_flag=True, help='Show what would be created without writing file')
@click.option('--no-validate', is_flag=True, help='Skip validation of sample existence')
@click.option('--output', '-o', type=click.Path(), help='Custom output path')
def create_assembly_collection(dataset_root, name, description, assembler, 
                               qc_filter, min_len, samples, sample_file,
                               co_assembly_collections, force, dry_run, 
                               no_validate, output):
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
        --qc-filter raw \\
        --samples sample1 --samples sample2 \\
        --co-assembly-collections tick_pools forest_pools
    
    \b
    # From file (samples.txt contains one sample per line)
    viral-snake create assembly-collection /path/to/dataset \\
        --name selected_samples \\
        --qc-filter raw__cutadapt_mgi_virome \\
        --sample-file samples.txt
    
    \b
    # Dry run to preview before creating
    viral-snake create assembly-collection /path/to/dataset \\
        --name test_collection \\
        --qc-filter raw \\
        --samples ALL \\
        --dry-run
    """
    
    dataset_path = Path(dataset_root)
    
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
            
            # Validate samples exist (unless skipped)
            if not no_validate:
                reads_dir = dataset_path / 'reads' / qc_filter
                valid_samples, missing_samples = validate_samples_exist(reads_dir, sample_list)
                
                if missing_samples:
                    console.print(f"[yellow]⚠ Warning:[/yellow] {len(missing_samples)} samples not found:")
                    for s in missing_samples[:5]:  # Show first 5
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
                'samples': sorted(sample_list)  # Sort for consistency
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
        'description': description or f'Assembly collection: {name}',
        'created': datetime.now().isoformat(),
        'created_by': f'viral-snake v{__version__}',
        'sources': sources,
        'stats': stats
    }
    
    # ========================================================================
    # Show preview
    # ========================================================================
    console.print("[bold]Preview:[/bold]\n")
    display_yaml_preview(collection, title=f"Assembly Collection: {name}")
    console.print()
    
    # Determine output path
    if output:
        config_file = Path(output)
    else:
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
    config_file.parent.mkdir(parents=True, exist_ok=True)
    with open(config_file, 'w') as f:
        yaml.dump(collection, f, default_flow_style=False, sort_keys=False)
    
    # ========================================================================
    # Success summary
    # ========================================================================
    console.print(Panel.fit(
        f"[bold green]✓ Assembly Collection Created[/bold green]\n\n"
        f"[cyan]File:[/cyan] {config_file}\n"
        f"[cyan]Name:[/cyan] {name}\n"
        f"[cyan]Individual assemblies:[/cyan] {stats['total_individual_assemblies']}\n"
        f"[cyan]Co-assemblies:[/cyan] {stats['total_co_assemblies']}",
        border_style="green"
    ))
    
    # ========================================================================
    # Next steps
    # ========================================================================
    console.print("\n[bold]Next steps:[/bold]")
    console.print(f"  1. Validate: [cyan]viral-snake validate dataset {dataset_root}[/cyan]")
    console.print(f"  2. Run DIAMOND: [cyan]viral-snake run diamond-batch {name}[/cyan]")
    console.print()
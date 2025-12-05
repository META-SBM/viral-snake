"""Create command group and shared utilities"""

import click
from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from rich.prompt import Confirm
from rich.progress import Progress, SpinnerColumn, TextColumn
from pathlib import Path
import yaml
import re

console = Console()


@click.group(name='create')
def create_group():
    """📝 Create metadata and configuration files"""
    pass


# ============================================================================
# Shared Helper Functions
# ============================================================================

def validate_collection_name(name):
    """Validate collection name contains only safe characters
    
    Args:
        name: Collection name to validate
        
    Raises:
        click.ClickException: If name contains invalid characters
    """
    if not re.match(r'^[a-zA-Z0-9_-]+$', name):
        raise click.ClickException(
            f"Invalid collection name: '{name}'\n"
            "Name can only contain: letters, numbers, underscores, and hyphens"
        )


def display_yaml_preview(data, title="Preview"):
    """Display YAML data with syntax highlighting"""
    yaml_str = yaml.dump(data, default_flow_style=False, sort_keys=False)
    syntax = Syntax(yaml_str, "yaml", theme="monokai", line_numbers=True)
    console.print(Panel(syntax, title=f"📄 {title}", border_style="cyan"))


def resolve_samples(dataset_path, qc_filter, samples, sample_file):
    """Resolve sample list from various input methods
    
    Args:
        dataset_path: Path to dataset root
        qc_filter: QC filter name
        samples: Tuple of sample names (can include 'ALL')
        sample_file: Path to file with sample names
        
    Returns:
        list: Resolved sample names
        
    Raises:
        click.ClickException: If samples can't be resolved
    """
    from ...utils import get_samples
    
    reads_dir = dataset_path / 'reads' / qc_filter
    
    if not reads_dir.exists():
        available = [d.name for d in (dataset_path / 'reads').iterdir() if d.is_dir()]
        raise click.ClickException(
            f"QC filter directory not found: {reads_dir}\n"
            f"Available filters: {', '.join(available)}"
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
            sample_list = [
                line.strip() for line in f 
                if line.strip() and not line.startswith('#')
            ]
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


def apply_exclusions(sample_list, exclude, exclude_file):
    """Apply exclusions to sample list
    
    Args:
        sample_list: List of samples
        exclude: Tuple of samples to exclude
        exclude_file: Path to file with exclusions
        
    Returns:
        tuple: (filtered_samples, excluded_samples, not_found_exclusions)
    """
    initial_count = len(sample_list)
    exclude_set = set()
    
    # From --exclude options
    if exclude:
        exclude_set.update(exclude)
        console.print(f"[yellow]Excluding {len(exclude)} samples specified via --exclude[/yellow]")
    
    # From --exclude-file
    if exclude_file:
        console.print(f"[yellow]Reading exclusions from:[/yellow] {exclude_file}")
        with open(exclude_file) as f:
            file_exclusions = [
                line.strip() for line in f 
                if line.strip() and not line.startswith('#')
            ]
            exclude_set.update(file_exclusions)
            console.print(f"[yellow]Loaded {len(file_exclusions)} exclusions from file[/yellow]")
    
    if not exclude_set:
        return sample_list, [], []
    
    # Track which exclusions matched
    excluded_samples = [s for s in sample_list if s in exclude_set]
    not_found = [e for e in exclude_set if e not in sample_list]
    
    # Apply filter
    filtered = [s for s in sample_list if s not in exclude_set]
    
    # Report
    console.print(f"\n[bold]Exclusion Summary:[/bold]")
    console.print(f"  Before exclusion: [cyan]{initial_count}[/cyan]")
    console.print(f"  Excluded: [yellow]{len(excluded_samples)}[/yellow]")
    console.print(f"  After exclusion: [green]{len(filtered)}[/green]")
    
    if excluded_samples:
        console.print(f"\n[dim]Excluded samples (showing first 10):[/dim]")
        for i, s in enumerate(sorted(excluded_samples)[:10], 1):
            console.print(f"  [dim]{i:3d}.[/dim] [yellow]{s}[/yellow]")
        if len(excluded_samples) > 10:
            console.print(f"  [dim]... and {len(excluded_samples) - 10} more[/dim]")
    
    if not_found:
        console.print(f"\n[yellow]⚠ Warning:[/yellow] {len(not_found)} exclusions didn't match:")
        for s in sorted(not_found)[:5]:
            console.print(f"  [dim]• {s}[/dim]")
        if len(not_found) > 5:
            console.print(f"  [dim]... and {len(not_found) - 5} more[/dim]")
    
    console.print()
    
    return filtered, excluded_samples, not_found


def validate_samples_exist(reads_dir, sample_list, show_progress=True):
    """Validate that all samples have R1 files
    
    Args:
        reads_dir: Directory containing reads
        sample_list: List of sample names to validate
        show_progress: Whether to show progress bar
        
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


# Import and register subcommands
from .feature_table import create_feature_table
from .assembly_collection import create_assembly_collection
from .sample_collection import create_sample_collection

create_group.add_command(create_feature_table)
create_group.add_command(create_assembly_collection)
create_group.add_command(create_sample_collection)
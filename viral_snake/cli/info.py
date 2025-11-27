"""Information and listing commands"""

import click
from rich.console import Console
from rich.table import Table
from pathlib import Path

from ..utils import get_samples

console = Console()


@click.group(name='info')
def info_group():
    """📊 Information about datasets and samples"""
    pass


@info_group.command(name='samples')
@click.argument('dataset_root', type=click.Path(exists=True))
@click.option('--qc-filter', default='raw', help='QC filter to use (default: raw)')
@click.option('--r1-pattern', default='_R1.fastq.gz', help='R1 file pattern')
@click.option('--r2-pattern', default='_R2.fastq.gz', help='R2 file pattern')
def list_samples(dataset_root, qc_filter, r1_pattern, r2_pattern):
    """List all samples in dataset"""
    
    dataset_path = Path(dataset_root)
    reads_dir = dataset_path / "reads" / qc_filter
    
    if not reads_dir.exists():
        console.print(f"[bold red]✗ Error:[/bold red] Reads directory not found: {reads_dir}")
        console.print(f"\n[yellow]Available QC filters:[/yellow]")
        
        reads_base = dataset_path / "reads"
        if reads_base.exists():
            for qc_dir in sorted(reads_base.iterdir()):
                if qc_dir.is_dir():
                    console.print(f"  • {qc_dir.name}")
        
        raise click.Abort()
    
    console.print(f"\n[bold cyan]Dataset:[/bold cyan] {dataset_root}")
    console.print(f"[bold cyan]QC Filter:[/bold cyan] {qc_filter}")
    console.print(f"[dim]Reads directory: {reads_dir}[/dim]")
    console.print(f"[dim]R1 pattern: {r1_pattern}[/dim]")
    console.print(f"[dim]R2 pattern: {r2_pattern}[/dim]\n")
    
    try:
        samples = get_samples(reads_dir, r1_pattern, r2_pattern)
        
        table = Table(
            title=f"📋 Discovered Samples ({qc_filter})", 
            show_header=True, 
            header_style="bold magenta"
        )
        table.add_column("#", style="dim", width=6, justify="right")
        table.add_column("Sample ID", style="cyan")
        table.add_column("Status", justify="center")
        
        for idx, sample in enumerate(samples, 1):
            r1_file = reads_dir / f"{sample}{r1_pattern}"
            r2_file = reads_dir / f"{sample}{r2_pattern}"
            
            r1_size = r1_file.stat().st_size if r1_file.exists() else 0
            r2_size = r2_file.stat().st_size if r2_file.exists() else 0
            
            status = "[green]✓[/green]" if (r1_size > 0 and r2_size > 0) else "[red]✗[/red]"
            
            table.add_row(str(idx), sample, status)
        
        console.print(table)
        console.print(f"\n[bold green]✓[/bold green] Found {len(samples)} sample(s) in [cyan]{qc_filter}[/cyan]\n")
        
    except FileNotFoundError as e:
        console.print(f"[bold red]✗ Error:[/bold red] {e}")
        raise click.Abort()


@info_group.command(name='qc-filters')
@click.argument('dataset_root', type=click.Path(exists=True))
def list_qc_filters(dataset_root):
    """List all available QC filters"""
    
    dataset_path = Path(dataset_root)
    reads_base = dataset_path / "reads"
    
    if not reads_base.exists():
        console.print(f"[bold red]✗ Error:[/bold red] Reads directory not found: {reads_base}")
        raise click.Abort()
    
    console.print(f"\n[bold cyan]Dataset:[/bold cyan] {dataset_root}\n")
    
    qc_filters = []
    for qc_dir in sorted(reads_base.iterdir()):
        if qc_dir.is_dir():
            try:
                samples = get_samples(qc_dir)
                qc_filters.append((qc_dir.name, len(samples)))
            except:
                qc_filters.append((qc_dir.name, 0))
    
    table = Table(title="🔬 Available QC Filters", show_header=True, header_style="bold magenta")
    table.add_column("QC Filter", style="cyan")
    table.add_column("Samples", justify="right", style="green")
    table.add_column("Default", justify="center")
    
    for qc_name, sample_count in qc_filters:
        is_default = "⭐" if qc_name == "raw" else ""
        table.add_row(qc_name, str(sample_count), is_default)
    
    console.print(table)
    console.print(f"\n[dim]Use --qc-filter to select a specific filter[/dim]\n")


@info_group.command(name='collections')
@click.argument('dataset_root', type=click.Path(exists=True))
def list_collections(dataset_root):
    """List all sample collections"""
    # TODO: implement
    console.print("[yellow]TODO: Implement collection listing[/yellow]")


@info_group.command(name='assemblies')
@click.argument('dataset_root', type=click.Path(exists=True))
def list_assemblies(dataset_root):
    """List all assemblies (individual + co-assemblies)"""
    # TODO: implement
    console.print("[yellow]TODO: Implement assembly listing[/yellow]")
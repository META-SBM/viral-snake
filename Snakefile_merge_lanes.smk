import glob
import os
import sys
from collections import defaultdict
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
from rich.columns import Columns

# Run identifier
RUN_NAME = "V350375304"

# Base paths
BASE_DIR = f"/mnt/mgx/DATASETS/INTERNAL/VIROME/VIROME5/{RUN_NAME}"
L01_DIR = os.path.join(BASE_DIR, "L01")
L02_DIR = os.path.join(BASE_DIR, "L02")
L03_DIR = os.path.join(BASE_DIR, "L03")
L04_DIR = os.path.join(BASE_DIR, "L04")
MERGED_DIR = os.path.join(BASE_DIR, "MERGED_LANES")

# Auto-detect sample IDs from L01 files
L01_files = glob.glob(os.path.join(L01_DIR, f"{RUN_NAME}_L01_*_1.fq.gz"))
ALL_SAMPLES = [os.path.basename(f).replace(f"{RUN_NAME}_L01_", "").replace("_1.fq.gz", "") for f in L01_files]

# Validate all lane files exist
lanes = ['L01', 'L02', 'L03', 'L04']
lane_dirs = [L01_DIR, L02_DIR, L03_DIR, L04_DIR]
issues_by_sample = defaultdict(lambda: {'R1': {'present': [], 'missing': []}, 'R2': {'present': [], 'missing': []}})

def get_file_size(filepath):
    """Get human-readable file size"""
    size = os.path.getsize(filepath)
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} TB"

samples_with_missing = set()
for sample in ALL_SAMPLES:
    has_missing = False
    for read in [1, 2]:
        read_key = f'R{read}'
        for lane, lane_dir in zip(lanes, lane_dirs):
            filepath = os.path.join(lane_dir, f"{RUN_NAME}_{lane}_{sample}_{read}.fq.gz")
            if os.path.exists(filepath):
                issues_by_sample[sample][read_key]['present'].append((filepath, get_file_size(filepath)))
            else:
                issues_by_sample[sample][read_key]['missing'].append(filepath)
                has_missing = True
    
    if has_missing:
        samples_with_missing.add(sample)
    else:
        # Remove samples with no issues from the issues dict
        if sample in issues_by_sample:
            del issues_by_sample[sample]

# Filter to only include samples with all files present
SAMPLES = [s for s in ALL_SAMPLES if s not in samples_with_missing]

if issues_by_sample:
    console = Console(stderr=True)
    console.print()
    console.print(Panel("[bold yellow]Warning: Some samples have missing lane files and will be skipped![/bold yellow]", 
                        title="WARNING", 
                        border_style="yellow bold",
                        box=box.DOUBLE))
    console.print()
    
    total_missing = 0
    for sample in sorted(issues_by_sample.keys()):
        
        # Create tables for R1 and R2
        tables_output = []
        for read in ['R1', 'R2']:
            table = Table(show_header=True, box=box.ROUNDED, title_style="bold yellow")
            table.add_column("Status", style="bold", width=10)
            table.add_column("File", style="cyan")
            table.add_column("Size", justify="right", style="dim")
            
            # Show present files
            for filepath, size in issues_by_sample[sample][read]['present']:
                table.add_row("[green]✓ Present[/green]", filepath, size)
            
            # Show missing files
            for filepath in issues_by_sample[sample][read]['missing']:
                table.add_row("[red]✗ Missing[/red]", filepath, "—")
                total_missing += 1
            
            tables_output.append(Panel(table, title=f"[bold yellow]{read}[/bold yellow]", border_style="yellow"))
        
        # Combine both tables in a panel
        console.print(Panel(
            Columns(tables_output, equal=True, expand=True),
            title=f"[bold cyan]Sample: {sample} (SKIPPED)[/bold cyan]",
            border_style="cyan",
            box=box.HEAVY
        ))
        console.print()
    
    console.print(Panel(
        f"[bold yellow]Skipping {len(issues_by_sample)} samples with {total_missing} missing files\n"
        f"Processing {len(SAMPLES)} samples with complete data[/bold yellow]",
        border_style="yellow bold",
        box=box.DOUBLE_EDGE
    ))
    console.print()

if not SAMPLES:
    console = Console(stderr=True)
    console.print(Panel("[bold red]ERROR: No samples with complete data found![/bold red]", 
                        border_style="red bold",
                        box=box.DOUBLE))
    sys.exit(1)

rule all:
    input:
        expand(os.path.join(MERGED_DIR, f"{RUN_NAME}_{{sample}}_{{read}}.fq.gz"), 
               sample=SAMPLES, read=[1, 2])

rule merge_lanes:
    input:
        l01 = os.path.join(L01_DIR, f"{RUN_NAME}_L01_{{sample}}_{{read}}.fq.gz"),
        l02 = os.path.join(L02_DIR, f"{RUN_NAME}_L02_{{sample}}_{{read}}.fq.gz"),
        l03 = os.path.join(L03_DIR, f"{RUN_NAME}_L03_{{sample}}_{{read}}.fq.gz"),
        l04 = os.path.join(L04_DIR, f"{RUN_NAME}_L04_{{sample}}_{{read}}.fq.gz")
    output:
        os.path.join(MERGED_DIR, f"{RUN_NAME}_{{sample}}_{{read}}.fq.gz")
    shell:
        "cat {input.l01} {input.l02} {input.l03} {input.l04} > {output}"
import glob
import os
import sys
from collections import defaultdict
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
from rich.columns import Columns

# ============================================================================
# CONFIGURATION - EDIT THESE FOR EACH RUN
# ============================================================================

# Run identifier
RUN_NAME = "V350375230"

# Specify which lanes to merge (only these will be checked and merged)
# Examples:
#   ["L02", "L03"]           # Only lanes 2 and 3
#   ["L01", "L02"]           # Only lanes 1 and 2
#   ["L01", "L02", "L03", "L04"]  # All four lanes
LANES_TO_MERGE = ["L01", "L02", "L03", "L04"]

# Base paths
BASE_DIR = f"/mnt/mgx/DATASETS/INTERNAL/VIROME/VIROME8/{RUN_NAME}"
MERGED_DIR = os.path.join(BASE_DIR, "MERGED_LANES")

# Optional: Path to sample sheet (set to None to auto-detect from files)
# SAMPLE_SHEET = "/path/to/your/sample_sheet.csv"
SAMPLE_SHEET = None

# ============================================================================
# AUTO-DETECTION AND VALIDATION
# ============================================================================

# Build lane directory paths based on configuration
LANE_DIRS = {lane: os.path.join(BASE_DIR, lane) for lane in LANES_TO_MERGE}

# Verify lane directories exist
missing_lane_dirs = [lane for lane, path in LANE_DIRS.items() if not os.path.exists(path)]
if missing_lane_dirs:
    console = Console(stderr=True)
    console.print(Panel(
        f"[bold red]ERROR: Lane directories not found: {', '.join(missing_lane_dirs)}[/bold red]\n"
        f"Expected paths: {[LANE_DIRS[l] for l in missing_lane_dirs]}",
        title="Configuration Error",
        border_style="red bold",
        box=box.DOUBLE
    ))
    sys.exit(1)

# Auto-detect sample IDs from first lane directory
first_lane = LANES_TO_MERGE[0]
first_lane_dir = LANE_DIRS[first_lane]
pattern = os.path.join(first_lane_dir, f"{RUN_NAME}_{first_lane}_*_1.fq.gz")
first_lane_files = glob.glob(pattern)

if not first_lane_files:
    console = Console(stderr=True)
    console.print(Panel(
        f"[bold red]ERROR: No files found in {first_lane} directory[/bold red]\n"
        f"Pattern searched: {pattern}",
        border_style="red bold",
        box=box.DOUBLE
    ))
    sys.exit(1)

ALL_SAMPLES = [
    os.path.basename(f).replace(f"{RUN_NAME}_{first_lane}_", "").replace("_1.fq.gz", "")
    for f in first_lane_files
]

# Optional: Filter by sample sheet
if SAMPLE_SHEET and os.path.exists(SAMPLE_SHEET):
    import pandas as pd
    console = Console(stderr=True)
    console.print(f"[cyan]Reading sample sheet: {SAMPLE_SHEET}[/cyan]")
    df = pd.read_csv(SAMPLE_SHEET)
    # Assuming 'ID' column contains sample names
    expected_samples = set(df['ID'].dropna().unique())
    console.print(f"[cyan]Expected {len(expected_samples)} samples from sheet[/cyan]")
    
    # Check for samples in files but not in sheet
    unexpected = set(ALL_SAMPLES) - expected_samples
    if unexpected:
        console.print(f"[yellow]Warning: {len(unexpected)} samples in files but not in sheet: {unexpected}[/yellow]")
    
    # Filter to only expected samples
    ALL_SAMPLES = [s for s in ALL_SAMPLES if s in expected_samples]

# Validate all required lane files exist
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
        for lane in LANES_TO_MERGE:
            lane_dir = LANE_DIRS[lane]
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

# ============================================================================
# REPORTING
# ============================================================================

console = Console(stderr=True)
console.print()
console.print(Panel(
    f"[bold cyan]Run: {RUN_NAME}[/bold cyan]\n"
    f"[bold cyan]Merging lanes: {', '.join(LANES_TO_MERGE)}[/bold cyan]\n"
    f"[dim]Detected {len(ALL_SAMPLES)} samples from {first_lane}[/dim]",
    title="Configuration",
    border_style="cyan",
    box=box.ROUNDED
))
console.print()

if issues_by_sample:
    console.print(Panel(
        "[bold yellow]Warning: Some samples have missing lane files and will be skipped![/bold yellow]", 
        title="WARNING", 
        border_style="yellow bold",
        box=box.DOUBLE
    ))
    console.print()
    
    total_missing = 0
    for sample in sorted(issues_by_sample.keys()):
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
else:
    console.print(Panel(
        f"[bold green]✓ All {len(SAMPLES)} samples have complete data in all {len(LANES_TO_MERGE)} lanes[/bold green]",
        border_style="green",
        box=box.ROUNDED
    ))
    console.print()

if not SAMPLES:
    console.print(Panel(
        "[bold red]ERROR: No samples with complete data found![/bold red]", 
        border_style="red bold",
        box=box.DOUBLE
    ))
    sys.exit(1)

# ============================================================================
# SNAKEMAKE RULES
# ============================================================================

rule all:
    input:
        expand(os.path.join(MERGED_DIR, f"{RUN_NAME}_{{sample}}_{{read}}.fq.gz"), 
               sample=SAMPLES, read=[1, 2])

def get_lane_inputs(wildcards):
    """Dynamically generate input files based on configured lanes"""
    return [
        os.path.join(LANE_DIRS[lane], f"{RUN_NAME}_{lane}_{wildcards.sample}_{wildcards.read}.fq.gz")
        for lane in LANES_TO_MERGE
    ]

rule merge_lanes:
    input:
        get_lane_inputs
    output:
        os.path.join(MERGED_DIR, f"{RUN_NAME}_{{sample}}_{{read}}.fq.gz")
    shell:
        "cat {input} > {output}"

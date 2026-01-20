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
# Filename Format:
# ============================================================================
# V350282721_SampleA__L01_L02_L03_L04_1.fq.gz
#          │        │  └──────────────┘ │
#          │        │         │         └─ Read (1 or 2)
#          │        │         └─────────── Technical: Lanes merged
#          │        └───────────────────── Double underscore separator
#          └────────────────────────────── Run name + Sample ID

# ============================================================================
# CONFIGURATION - EDIT THESE FOR EACH RUN
# ============================================================================

# File system and dataset configuration
FS_PREFIX = "/mnt/mgx/DATASETS/INTERNAL/VIROME/"
DATASET = "VIROME2"

# Run identifier
RUN_NAME = "V350282721"

# Specify which lanes to merge (only these will be checked and merged)
# Examples:
#   ["L02", "L03"]           # Only lanes 2 and 3
#   ["L01", "L02"]           # Only lanes 1 and 2
#   ["L01", "L02", "L03", "L04"]  # All four lanes
LANES_TO_MERGE = ["L01", "L02", "L03", "L04"]
LANES_TO_MERGE = ["L01", "L02"]
LANES_TO_MERGE = ["L03", "L04"]

# Create lane string for filenames (e.g., "L01_L02_L03_L04")
LANE_STRING = "_".join(LANES_TO_MERGE)

# Build paths from configuration
BASE_DIR = os.path.join(FS_PREFIX, DATASET, RUN_NAME)
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
    f"[bold cyan]File System: {FS_PREFIX}[/bold cyan]\n"
    f"[bold cyan]Dataset: {DATASET}[/bold cyan]\n"
    f"[bold cyan]Run: {RUN_NAME}[/bold cyan]\n"
    f"[bold cyan]Merging lanes: {', '.join(LANES_TO_MERGE)}[/bold cyan]\n"
    f"[bold cyan]Output naming: {RUN_NAME}_{{sample}}__{LANE_STRING}_{{read}}.fq.gz[/bold cyan]\n"
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
# STATS REPORT GENERATION
# ============================================================================

def generate_stats_report():
    """Generate a comprehensive stats report with improved formatting"""
    import json
    from datetime import datetime
    
    stats_file = os.path.join(BASE_DIR, f"{RUN_NAME}_merge_stats.txt")
    json_file = os.path.join(BASE_DIR, f"{RUN_NAME}_merge_stats.json")
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Prepare structured data for JSON
    stats_data = {
        "timestamp": timestamp,
        "filesystem_prefix": FS_PREFIX,
        "dataset": DATASET,
        "run_name": RUN_NAME,
        "lanes_merged": LANES_TO_MERGE,
        "lane_string": LANE_STRING,
        "base_directory": BASE_DIR,
        "output_directory": MERGED_DIR,
        "output_naming_pattern": f"{RUN_NAME}_{{sample}}__{LANE_STRING}_{{read}}.fq.gz",
        "sample_sheet": SAMPLE_SHEET,
        "total_samples_detected": len(ALL_SAMPLES),
        "samples_to_process": len(SAMPLES),
        "samples_skipped": len(issues_by_sample),
        "processed_samples": sorted(SAMPLES),
        "skipped_samples": {}
    }
    
    # Build detailed skip information
    for sample in sorted(issues_by_sample.keys()):
        stats_data["skipped_samples"][sample] = {
            "R1_present": [f[0] for f in issues_by_sample[sample]['R1']['present']],
            "R1_missing": issues_by_sample[sample]['R1']['missing'],
            "R2_present": [f[0] for f in issues_by_sample[sample]['R2']['present']],
            "R2_missing": issues_by_sample[sample]['R2']['missing']
        }
    
    # Write human-readable text report with improved formatting
    with open(stats_file, 'w') as f:
        # Header
        f.write("\n")
        f.write("╔" + "═" * 78 + "╗\n")
        f.write("║" + " " * 78 + "║\n")
        f.write("║" + "       LANE MERGE STATISTICS REPORT".center(78) + "║\n")
        f.write("║" + " " * 78 + "║\n")
        f.write("╚" + "═" * 78 + "╝\n")
        f.write("\n")
        
        # Configuration Box
        f.write("┌─ CONFIGURATION " + "─" * 61 + "┐\n")
        f.write("│                                                                              │\n")
        f.write(f"│  Generated      : {timestamp:<59} │\n")
        f.write(f"│  FS Prefix      : {FS_PREFIX:<59} │\n")
        f.write(f"│  Dataset        : {DATASET:<59} │\n")
        f.write(f"│  Run Name       : {RUN_NAME:<59} │\n")
        f.write(f"│  Lanes Merged   : {', '.join(LANES_TO_MERGE):<59} │\n")
        f.write(f"│  Lane String    : {LANE_STRING:<59} │\n")
        f.write(f"│  Base Directory : {BASE_DIR:<59} │\n")
        f.write(f"│  Output Dir     : {MERGED_DIR:<59} │\n")
        output_pattern = f"{RUN_NAME}_{{sample}}__{LANE_STRING}_{{read}}.fq.gz"
        f.write(f"│  Output Pattern : {output_pattern:<59} │\n")
        sample_sheet_display = SAMPLE_SHEET if SAMPLE_SHEET else "None (auto-detected)"
        f.write(f"│  Sample Sheet   : {sample_sheet_display:<59} │\n")
        f.write("│                                                                              │\n")
        f.write("└" + "─" * 78 + "┘\n")
        f.write("\n\n")
        
        # Summary Box
        f.write("┌─ SUMMARY " + "─" * 67 + "┐\n")
        f.write("│                                                                              │\n")
        f.write(f"│  Total Samples Detected  : {len(ALL_SAMPLES):>3}                                                   │\n")
        f.write(f"│  Samples to Process      : {len(SAMPLES):>3}  ✓                                                │\n")
        f.write(f"│  Samples Skipped         : {len(issues_by_sample):>3}  ✗                                                │\n")
        f.write("│                                                                              │\n")
        f.write("└" + "─" * 78 + "┘\n")
        f.write("\n\n")
        
        # Samples to Process
        if SAMPLES:
            f.write("╔" + "═" * 78 + "╗\n")
            f.write(f"║  SAMPLES TO PROCESS ({len(SAMPLES)})".ljust(79) + "║\n")
            f.write("╚" + "═" * 78 + "╝\n")
            f.write("\n")
            
            for idx, sample in enumerate(sorted(SAMPLES), 1):
                f.write(f"[{idx:>3}] {sample}\n")
                f.write("     " + "─" * 70 + "\n")
                
                # R1 and R2 in columns
                f.write("     READ 1 (R1)                              READ 2 (R2)\n")
                f.write("     " + "─" * 70 + "\n")
                
                # Input files
                f.write("     Inputs:\n")
                for lane_idx, lane in enumerate(LANES_TO_MERGE):
                    r1_file = os.path.join(LANE_DIRS[lane], f"{RUN_NAME}_{lane}_{sample}_1.fq.gz")
                    r2_file = os.path.join(LANE_DIRS[lane], f"{RUN_NAME}_{lane}_{sample}_2.fq.gz")
                    
                    # Get file sizes
                    r1_size = get_file_size(r1_file)
                    r2_size = get_file_size(r2_file)
                    
                    f.write(f"       [{lane}] {r1_size:>8}                         [{lane}] {r2_size:>8}\n")
                
                # Output files with lane info in filename
                f.write("\n")
                f.write("     Outputs:\n")
                r1_output = os.path.join(MERGED_DIR, f"{RUN_NAME}_{sample}__{LANE_STRING}_1.fq.gz")
                r2_output = os.path.join(MERGED_DIR, f"{RUN_NAME}_{sample}__{LANE_STRING}_2.fq.gz")
                f.write(f"       → {os.path.basename(r1_output):<30} → {os.path.basename(r2_output)}\n")
                
                f.write("\n")
        
        # Samples Skipped
        if issues_by_sample:
            f.write("\n")
            f.write("╔" + "═" * 78 + "╗\n")
            f.write(f"║  SAMPLES SKIPPED ({len(issues_by_sample)})".ljust(79) + "║\n")
            f.write("╚" + "═" * 78 + "╝\n")
            f.write("\n")
            
            for idx, sample in enumerate(sorted(issues_by_sample.keys()), 1):
                f.write(f"[{idx:>3}] {sample}  ✗ SKIPPED\n")
                f.write("     " + "─" * 70 + "\n")
                f.write("     Reason: Missing lane files\n")
                f.write("\n")
                
                # Count missing files
                r1_missing_count = len(issues_by_sample[sample]['R1']['missing'])
                r2_missing_count = len(issues_by_sample[sample]['R2']['missing'])
                r1_present_count = len(issues_by_sample[sample]['R1']['present'])
                r2_present_count = len(issues_by_sample[sample]['R2']['present'])
                
                f.write(f"     READ 1: {r1_present_count}/{len(LANES_TO_MERGE)} lanes present, {r1_missing_count} missing\n")
                f.write(f"     READ 2: {r2_present_count}/{len(LANES_TO_MERGE)} lanes present, {r2_missing_count} missing\n")
                f.write("\n")
                
                # Show detailed status for each lane
                f.write("     Lane Status:\n")
                f.write("     ┌" + "─" * 68 + "┐\n")
                f.write("     │ Lane  │ R1 Status          │ R2 Status          │ Details     │\n")
                f.write("     ├" + "─" * 68 + "┤\n")
                
                for lane in LANES_TO_MERGE:
                    # Check R1 status
                    r1_path = os.path.join(LANE_DIRS[lane], f"{RUN_NAME}_{lane}_{sample}_1.fq.gz")
                    r2_path = os.path.join(LANE_DIRS[lane], f"{RUN_NAME}_{lane}_{sample}_2.fq.gz")
                    
                    r1_present = any(r1_path == f[0] for f in issues_by_sample[sample]['R1']['present'])
                    r2_present = any(r2_path == f[0] for f in issues_by_sample[sample]['R2']['present'])
                    
                    r1_status = "✓ Present" if r1_present else "✗ MISSING"
                    r2_status = "✓ Present" if r2_present else "✗ MISSING"
                    
                    # Get sizes if present
                    r1_detail = ""
                    r2_detail = ""
                    if r1_present:
                        r1_detail = next((f[1] for f in issues_by_sample[sample]['R1']['present'] if f[0] == r1_path), "")
                    if r2_present:
                        r2_detail = next((f[1] for f in issues_by_sample[sample]['R2']['present'] if f[0] == r2_path), "")
                    
                    detail = f"{r1_detail}/{r2_detail}" if (r1_detail or r2_detail) else "—"
                    
                    f.write(f"     │  {lane}  │ {r1_status:<18} │ {r2_status:<18} │ {detail:<11} │\n")
                
                f.write("     └" + "─" * 68 + "┘\n")
                f.write("\n")
                
                # Show missing file paths
                all_missing = issues_by_sample[sample]['R1']['missing'] + issues_by_sample[sample]['R2']['missing']
                if all_missing:
                    f.write("     Missing Files:\n")
                    for missing_path in all_missing:
                        f.write(f"       ✗ {missing_path}\n")
                    f.write("\n")
        
        # Footer
        f.write("\n")
        f.write("─" * 80 + "\n")
        f.write(f"End of Report | Generated: {timestamp}\n")
        f.write("─" * 80 + "\n")
    
    # Write JSON report
    with open(json_file, 'w') as f:
        json.dump(stats_data, f, indent=2)
    
    console.print(Panel(
        f"[bold green]✓ Stats reports saved:[/bold green]\n"
        f"  • Text: {stats_file}\n"
        f"  • JSON: {json_file}",
        title="Reports Generated",
        border_style="green",
        box=box.ROUNDED
    ))

# Generate the report
generate_stats_report()

# ============================================================================
# SNAKEMAKE RULES
# ============================================================================

rule all:
    input:
        expand(os.path.join(MERGED_DIR, f"{RUN_NAME}_{{sample}}__{LANE_STRING}_{{read}}.fq.gz"), 
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
        os.path.join(MERGED_DIR, f"{RUN_NAME}_{{sample}}__{LANE_STRING}_{{read}}.fq.gz")
    shell:
        "cat {input} > {output}"
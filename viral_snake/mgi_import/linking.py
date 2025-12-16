import os
import glob
import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
from rich.progress import Progress, SpinnerColumn, TextColumn
from collections import defaultdict

def create_symlinks_multi_lane(
    df_renamed,
    run_name,
    base_dir,
    target_dir,
    merged_lanes=[],
    standalone_lanes=[],
    on_missing='warn',
    dry_run=True
):
    """
    Create symlinks from sequencing files to target directory with sample names.
    
    Parameters:
    -----------
    df_renamed : pd.DataFrame
        DataFrame with renamed columns (Barcode_L01, ID_L01, etc.)
    run_name : str
        Run identifier (e.g., "V350375302")
    base_dir : str
        Base directory containing L01, L02, etc. and MERGED_LANES subdirs
    target_dir : str
        Target directory where symlinks will be created
    merged_lanes : list
        List of lanes that were merged together (e.g., ['L02', 'L03'])
    standalone_lanes : list
        List of lanes that remain separate (e.g., ['L01', 'L04'])
    on_missing : str
        How to handle missing files: 'warn', 'error', or 'skip'
    dry_run : bool
        If True, only show what would be done without creating symlinks
    
    Returns:
    --------
    dict : Summary statistics of the operation
    """
    
    console = Console()
    
    # Header
    mode = "[yellow bold]DRY RUN MODE[/yellow bold]" if dry_run else "[green bold]EXECUTION MODE[/green bold]"
    console.print(Panel(
        f"{mode}\n\n"
        f"[cyan]Run:[/cyan] {run_name}\n"
        f"[cyan]Source:[/cyan] {base_dir}\n"
        f"[cyan]Target:[/cyan] {target_dir}\n"
        f"[cyan]Merged lanes:[/cyan] {', '.join(merged_lanes) if merged_lanes else 'None'}\n"
        f"[cyan]Standalone lanes:[/cyan] {', '.join(standalone_lanes) if standalone_lanes else 'None'}",
        title="🔗 Symlink Creator",
        box=box.DOUBLE,
        border_style="cyan"
    ))
    
    # Step 1: Verify merged lanes have identical mappings
    if len(merged_lanes) > 1:
        console.print(Panel(
            f"[cyan]Verifying {len(merged_lanes)} merged lanes have identical barcode→ID mappings...[/cyan]",
            box=box.ROUNDED
        ))
        
        reference_lane = merged_lanes[0]
        ref_barcode_col = f'Barcode_{reference_lane}'
        ref_id_col = f'ID_{reference_lane}'
        
        if ref_barcode_col not in df_renamed.columns or ref_id_col not in df_renamed.columns:
            console.print(f"[red]ERROR: Columns {ref_barcode_col} or {ref_id_col} not found![/red]")
            return None
        
        # Build reference mapping
        ref_mapping = df_renamed[[ref_barcode_col, ref_id_col]].dropna()
        ref_mapping = dict(zip(ref_mapping[ref_barcode_col], ref_mapping[ref_id_col]))
        
        # Compare with other merged lanes
        all_match = True
        for lane in merged_lanes[1:]:
            barcode_col = f'Barcode_{lane}'
            id_col = f'ID_{lane}'
            
            if barcode_col not in df_renamed.columns or id_col not in df_renamed.columns:
                console.print(f"[red]ERROR: Columns {barcode_col} or {id_col} not found![/red]")
                return None
            
            lane_mapping = df_renamed[[barcode_col, id_col]].dropna()
            lane_mapping = dict(zip(lane_mapping[barcode_col], lane_mapping[id_col]))
            
            # Check for differences
            differences = []
            all_barcodes = set(ref_mapping.keys()) | set(lane_mapping.keys())
            
            for barcode in all_barcodes:
                ref_id = ref_mapping.get(barcode)
                lane_id = lane_mapping.get(barcode)
                if ref_id != lane_id:
                    differences.append((barcode, ref_id, lane_id))
            
            if differences:
                all_match = False
                console.print(f"[red]✗ Mapping mismatch between {reference_lane} and {lane}![/red]")
                error_table = Table(title="Mapping Differences", box=box.ROUNDED, border_style="red")
                error_table.add_column("Barcode", style="cyan")
                error_table.add_column(reference_lane, style="yellow")
                error_table.add_column(lane, style="yellow")
                
                for barcode, ref_id, lane_id in differences[:10]:  # Show first 10
                    error_table.add_row(str(barcode), str(ref_id), str(lane_id))
                
                console.print(error_table)
                if len(differences) > 10:
                    console.print(f"[dim]... and {len(differences) - 10} more differences[/dim]")
                return None
        
        if all_match:
            console.print(f"[green]✓ All {len(merged_lanes)} merged lanes have identical mappings ({len(ref_mapping)} barcodes)[/green]\n")
    
    # Step 2: Build mappings for each source
    console.print("[cyan]Building barcode → sample ID mappings...[/cyan]")
    mappings = {}
    
    # Merged lanes mapping (use first lane as reference)
    if merged_lanes:
        ref_lane = merged_lanes[0]
        barcode_col = f'Barcode_{ref_lane}'
        id_col = f'ID_{ref_lane}'
        merged_df = df_renamed[[barcode_col, id_col]].dropna()
        mappings['MERGED'] = dict(zip(merged_df[barcode_col], merged_df[id_col]))
        console.print(f"  [green]✓[/green] MERGED ({', '.join(merged_lanes)}): {len(mappings['MERGED'])} samples")
    
    # Standalone lanes mappings
    for lane in standalone_lanes:
        barcode_col = f'Barcode_{lane}'
        id_col = f'ID_{lane}'
        
        if barcode_col not in df_renamed.columns or id_col not in df_renamed.columns:
            console.print(f"  [yellow]⚠[/yellow] {lane}: Columns not found, skipping")
            continue
        
        lane_df = df_renamed[[barcode_col, id_col]].dropna()
        mappings[lane] = dict(zip(lane_df[barcode_col], lane_df[id_col]))
        console.print(f"  [green]✓[/green] {lane}: {len(mappings[lane])} samples")
    
    console.print()
    
    # Step 3: Check file availability and validate pairs
    console.print("[cyan]Scanning source files...[/cyan]")
    
    file_status = {}
    
    for source, barcode_to_id in mappings.items():
        file_status[source] = {
            'valid': {},      # sample_id: {'R1': path, 'R2': path}
            'incomplete': {}, # sample_id: {'R1': path} or {'R2': path}
            'missing': []     # [(sample_id, barcode)]
        }
        
        # Determine source directory
        if source == 'MERGED':
            source_dir = os.path.join(base_dir, "MERGED_LANES")
        else:
            source_dir = os.path.join(base_dir, source)
        
        if not os.path.exists(source_dir):
            console.print(f"  [red]✗[/red] {source}: Directory not found: {source_dir}")
            continue
        
        # Check each sample
        for barcode, sample_id in barcode_to_id.items():
            reads_found = {}
            
            for read in [1, 2]:
                if source == 'MERGED':
                    source_file = os.path.join(source_dir, f"{run_name}_{barcode}_{read}.fq.gz")
                else:
                    source_file = os.path.join(source_dir, f"{run_name}_{source}_{barcode}_{read}.fq.gz")
                
                if os.path.exists(source_file):
                    reads_found[f'R{read}'] = source_file
            
            if len(reads_found) == 2:
                file_status[source]['valid'][sample_id] = reads_found
            elif len(reads_found) == 1:
                file_status[source]['incomplete'][sample_id] = reads_found
            else:
                file_status[source]['missing'].append((sample_id, barcode))
        
        # Report for this source
        n_valid = len(file_status[source]['valid'])
        n_incomplete = len(file_status[source]['incomplete'])
        n_missing = len(file_status[source]['missing'])
        
        status_msg = f"  {source}: "
        if n_valid > 0:
            status_msg += f"[green]{n_valid} complete pairs[/green]"
        if n_incomplete > 0:
            status_msg += f"  [yellow]{n_incomplete} incomplete[/yellow]"
        if n_missing > 0:
            status_msg += f"  [red]{n_missing} missing[/red]"
        
        console.print(status_msg)
    
    console.print()
    
    # Step 4: Show warnings for incomplete/missing samples
    total_incomplete = sum(len(file_status[s]['incomplete']) for s in file_status)
    total_missing = sum(len(file_status[s]['missing']) for s in file_status)
    
    if total_incomplete > 0:
        console.print(Panel(
            "[yellow]Warning: Some samples have incomplete read pairs (only R1 or R2 found)[/yellow]",
            title="⚠ Incomplete Samples",
            border_style="yellow"
        ))
        
        for source in file_status:
            if file_status[source]['incomplete']:
                console.print(f"[yellow]{source}:[/yellow]")
                for sample_id, reads in list(file_status[source]['incomplete'].items())[:5]:
                    missing_read = 'R2' if 'R1' in reads else 'R1'
                    console.print(f"  • {sample_id} (missing {missing_read})")
                
                if len(file_status[source]['incomplete']) > 5:
                    console.print(f"  [dim]... and {len(file_status[source]['incomplete']) - 5} more[/dim]")
        console.print()
    
    if total_missing > 0:
        if on_missing == 'error':
            console.print(Panel(
                "[red]ERROR: Missing files detected and on_missing='error'[/red]",
                title="❌ Error",
                border_style="red"
            ))
            return None
        elif on_missing == 'warn':
            console.print(Panel(
                "[yellow]Warning: Some samples have no files found (both R1 and R2 missing)[/yellow]",
                title="⚠ Missing Files",
                border_style="yellow"
            ))
            
            for source in file_status:
                if file_status[source]['missing']:
                    console.print(f"[yellow]{source}:[/yellow]")
                    for sample_id, barcode in list(file_status[source]['missing'])[:5]:
                        console.print(f"  • {sample_id} (barcode: {barcode})")
                    
                    if len(file_status[source]['missing']) > 5:
                        console.print(f"  [dim]... and {len(file_status[source]['missing']) - 5} more[/dim]")
            console.print()
    
    # Step 5: Create symlinks for valid samples
    if not dry_run:
        os.makedirs(target_dir, exist_ok=True)
    
    stats = {
        'total_samples': 0,
        'total_links': 0,
        'created': 0,
        'failed': 0,
        'by_source': {}
    }
    
    for source in file_status:
        valid_samples = file_status[source]['valid']
        
        if not valid_samples:
            continue
        
        # Create table for this source
        action = "Will create" if dry_run else "Creating"
        source_display = f"MERGED ({', '.join(merged_lanes)})" if source == 'MERGED' else source
        
        table = Table(
            title=f"{action} symlinks for {source_display}",
            box=box.ROUNDED,
            show_lines=False
        )
        table.add_column("Sample ID", style="cyan", no_wrap=True)
        table.add_column("R1", justify="center", width=12)
        table.add_column("R2", justify="center", width=12)
        
        source_stats = {'samples': len(valid_samples), 'created': 0, 'failed': 0}
        
        for sample_id in sorted(valid_samples.keys()):
            reads = valid_samples[sample_id]
            statuses = []
            
            for read in ['R1', 'R2']:
                source_file = reads[read]
                target_file = os.path.join(target_dir, f"{sample_id}_{read}.fastq.gz")
                
                if dry_run:
                    statuses.append("[yellow]pending[/yellow]")
                else:
                    try:
                        if os.path.islink(target_file) or os.path.exists(target_file):
                            os.remove(target_file)
                        os.symlink(source_file, target_file)
                        statuses.append("[green]✓[/green]")
                        source_stats['created'] += 1
                        stats['created'] += 1
                    except Exception as e:
                        statuses.append("[red]✗[/red]")
                        source_stats['failed'] += 1
                        stats['failed'] += 1
                        console.print(f"[red]Error: {sample_id} {read}: {e}[/red]")
            
            table.add_row(sample_id, statuses[0], statuses[1])
        
        console.print(table)
        console.print()
        
        stats['by_source'][source] = source_stats
        stats['total_samples'] += source_stats['samples']
        stats['total_links'] += source_stats['samples'] * 2
    
    # Step 6: Final summary
    summary_lines = []
    
    if dry_run:
        summary_lines.append(f"[yellow]Would create {stats['total_links']} symlinks for {stats['total_samples']} samples[/yellow]")
        summary_lines.append(f"[dim]Set dry_run=False to execute[/dim]")
    else:
        summary_lines.append(f"[green]Created: {stats['created']}[/green]  [red]Failed: {stats['failed']}[/red]")
        summary_lines.append(f"[dim]Processed {stats['total_samples']} samples across {len(stats['by_source'])} source(s)[/dim]")
    
    if total_incomplete > 0:
        summary_lines.append(f"[yellow]⚠ {total_incomplete} samples with incomplete pairs (skipped)[/yellow]")
    
    if total_missing > 0 and on_missing != 'skip':
        summary_lines.append(f"[yellow]⚠ {total_missing} samples with missing files[/yellow]")
    
    console.print(Panel(
        "\n".join(summary_lines),
        title="📊 Summary",
        box=box.DOUBLE,
        border_style="green" if not dry_run and stats['failed'] == 0 else "yellow"
    ))
    
    return stats
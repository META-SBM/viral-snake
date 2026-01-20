import os
import glob
import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
from rich.progress import Progress, SpinnerColumn, TextColumn
from collections import defaultdict

def get_file_size(filepath):
    """Get human-readable file size"""
    if not os.path.exists(filepath):
        return None
    size = os.path.getsize(filepath)
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} TB"

def get_file_size_bytes(filepath):
    """Get file size in bytes"""
    if not os.path.exists(filepath):
        return None
    return os.path.getsize(filepath)

def create_symlinks_multi_lane(
    df_renamed,
    run_name,
    base_dir,
    target_dir,
    merged_lanes=[],
    standalone_lanes=[],
    on_missing='warn',
    dry_run=True,
    min_file_size_mb=10
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
    min_file_size_mb : int
        Minimum file size in MB to flag as suspicious (default: 10)
    
    Returns:
    --------
    dict : Summary statistics of the operation
    """
    
    console = Console()
    
    # ═══════════════════════════════════════════════════════════════
    # HEADER
    # ═══════════════════════════════════════════════════════════════
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
    
    # ═══════════════════════════════════════════════════════════════
    # STEP 1: VERIFY MERGED LANES HAVE IDENTICAL MAPPINGS
    # ═══════════════════════════════════════════════════════════════
    if len(merged_lanes) > 1:
        console.print(Panel(
            f"[cyan]Verifying {len(merged_lanes)} merged lanes have identical barcode→ID mappings...[/cyan]",
            box=box.ROUNDED
        ))
        
        reference_lane = merged_lanes[0]
        ref_barcode_col = f'Barcode_{reference_lane}'
        ref_id_col = f'ID_{reference_lane}'
        
        if ref_barcode_col not in df_renamed.columns or ref_id_col not in df_renamed.columns:
            console.print(Panel(
                f"[red bold]ERROR: Required columns not found![/red bold]\n\n"
                f"Looking for: [yellow]{ref_barcode_col}[/yellow] and [yellow]{ref_id_col}[/yellow]\n"
                f"Available columns: [dim]{', '.join(df_renamed.columns)}[/dim]",
                title="❌ Column Error",
                border_style="red",
                box=box.DOUBLE
            ))
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
                console.print(Panel(
                    f"[red bold]ERROR: Required columns not found![/red bold]\n\n"
                    f"Looking for: [yellow]{barcode_col}[/yellow] and [yellow]{id_col}[/yellow]\n"
                    f"Available columns: [dim]{', '.join(df_renamed.columns)}[/dim]",
                    title="❌ Column Error",
                    border_style="red",
                    box=box.DOUBLE
                ))
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
                console.print(Panel(
                    f"[red bold]MAPPING MISMATCH![/red bold]\n\n"
                    f"Lanes [yellow]{reference_lane}[/yellow] and [yellow]{lane}[/yellow] have different barcode→ID mappings.\n"
                    f"Found [red]{len(differences)}[/red] mismatched entries.\n\n"
                    "[yellow]These lanes cannot be merged because they map barcodes to different sample IDs![/yellow]",
                    title="❌ Merge Validation Failed",
                    border_style="red",
                    box=box.DOUBLE
                ))
                
                error_table = Table(
                    title="Mapping Differences",
                    box=box.ROUNDED,
                    border_style="red",
                    show_lines=True
                )
                error_table.add_column("Barcode", style="cyan", no_wrap=True)
                error_table.add_column(reference_lane, style="yellow")
                error_table.add_column(lane, style="yellow")
                error_table.add_column("Issue", style="red dim")
                
                for barcode, ref_id, lane_id in differences[:10]:
                    if ref_id is None:
                        issue = "Missing in ref"
                    elif lane_id is None:
                        issue = "Missing in lane"
                    else:
                        issue = "Different IDs"
                    error_table.add_row(str(barcode), str(ref_id), str(lane_id), issue)
                
                console.print(error_table)
                if len(differences) > 10:
                    console.print(f"[dim]... and {len(differences) - 10} more differences[/dim]\n")
                
                return None
        
        if all_match:
            console.print(f"[green]✓ All {len(merged_lanes)} merged lanes have identical mappings ({len(ref_mapping)} barcodes)[/green]\n")
    
    # ═══════════════════════════════════════════════════════════════
    # STEP 2: BUILD MAPPINGS FOR EACH SOURCE
    # ═══════════════════════════════════════════════════════════════
    console.print("[cyan bold]Building barcode → sample ID mappings...[/cyan bold]")
    mappings = {}
    
    # Merged lanes mapping (use first lane as reference)
    if merged_lanes:
        ref_lane = merged_lanes[0]
        barcode_col = f'Barcode_{ref_lane}'
        id_col = f'ID_{ref_lane}'
        merged_df = df_renamed[[barcode_col, id_col]].dropna()
        mappings['MERGED'] = dict(zip(merged_df[barcode_col], merged_df[id_col]))
        console.print(f"  [green]✓[/green] MERGED ({', '.join(merged_lanes)}): [cyan]{len(mappings['MERGED'])}[/cyan] samples")
    
    # Standalone lanes mappings
    for lane in standalone_lanes:
        barcode_col = f'Barcode_{lane}'
        id_col = f'ID_{lane}'
        
        if barcode_col not in df_renamed.columns or id_col not in df_renamed.columns:
            console.print(f"  [yellow]⚠[/yellow] {lane}: Columns [yellow]{barcode_col}[/yellow] or [yellow]{id_col}[/yellow] not found, skipping")
            continue
        
        lane_df = df_renamed[[barcode_col, id_col]].dropna()
        mappings[lane] = dict(zip(lane_df[barcode_col], lane_df[id_col]))
        console.print(f"  [green]✓[/green] {lane}: [cyan]{len(mappings[lane])}[/cyan] samples")
    
    if not mappings:
        console.print(Panel(
            "[red bold]ERROR: No valid mappings found![/red bold]\n\n"
            "Check that your DataFrame has the correct column names:\n"
            "• Barcode_L01, ID_L01, etc. for standalone lanes\n"
            "• Columns must match the lanes specified in merged_lanes and standalone_lanes",
            title="❌ Configuration Error",
            border_style="red",
            box=box.DOUBLE
        ))
        return None
    
    console.print()
    
    # ═══════════════════════════════════════════════════════════════
    # STEP 3: CHECK FILE AVAILABILITY AND VALIDATE PAIRS
    # ═══════════════════════════════════════════════════════════════
    console.print("[cyan bold]Scanning source files...[/cyan bold]")
    
    file_status = {}
    min_file_size_bytes = min_file_size_mb * 1024 * 1024
    total_small_files = 0
    
    for source, barcode_to_id in mappings.items():
        file_status[source] = {
            'valid': {},      # sample_id: {'R1': path, 'R2': path}
            'incomplete': {}, # sample_id: {'found': {R1/R2: path}, 'lane_status': {...}}
            'missing': []     # [(sample_id, barcode, lane_status)]
        }
        
        # Determine which directory to check
        if source == 'MERGED':
            # For merged lanes, check the MERGED_LANES directory
            source_dir = os.path.join(base_dir, 'MERGED_LANES')
            lanes_to_check = ['MERGED_LANES']  # Single virtual lane
        else:
            # For standalone lanes, check the specific lane directory
            source_dir = os.path.join(base_dir, source)
            lanes_to_check = [source]
        
        # Check each sample
        for barcode, sample_id in barcode_to_id.items():
            # Build detailed lane status
            lane_status = {}
            
            for lane in lanes_to_check:
                if source == 'MERGED':
                    lane_dir = source_dir
                else:
                    lane_dir = os.path.join(base_dir, lane)
                
                if not os.path.exists(lane_dir):
                    lane_status[lane] = {
                        'R1': {'exists': False, 'path': None, 'size': None, 'size_bytes': None},
                        'R2': {'exists': False, 'path': None, 'size': None, 'size_bytes': None},
                        'dir_exists': False
                    }
                    continue
                
                lane_status[lane] = {'dir_exists': True}
                
                for read in ['R1', 'R2']:
                    read_num = read[1]  # '1' or '2'
                    if source == 'MERGED':
                        # For merged lanes, files are in MERGED_LANES with special naming
                        filepath = os.path.join(lane_dir, f"{run_name}_{barcode}_{read_num}.fq.gz")
                    else:
                        filepath = os.path.join(lane_dir, f"{run_name}_{lane}_{barcode}_{read_num}.fq.gz")
                    
                    if os.path.exists(filepath):
                        file_size = get_file_size(filepath)
                        file_size_bytes = get_file_size_bytes(filepath)
                        
                        lane_status[lane][read] = {
                            'exists': True,
                            'path': filepath,
                            'size': file_size,
                            'size_bytes': file_size_bytes
                        }
                        
                        # Check for suspiciously small files
                        if file_size_bytes and file_size_bytes < min_file_size_bytes:
                            lane_status[lane][read]['is_small'] = True
                            total_small_files += 1
                        else:
                            lane_status[lane][read]['is_small'] = False
                    else:
                        lane_status[lane][read] = {
                            'exists': False,
                            'path': filepath,
                            'size': None,
                            'size_bytes': None,
                            'is_small': False
                        }
            
            # FIX: Check if any lanes exist before using all()
            existing_lanes = [lane for lane in lanes_to_check if lane_status[lane].get('dir_exists', False)]
            
            if not existing_lanes:
                # All directories missing - sample is completely missing
                file_status[source]['missing'].append((sample_id, barcode, lane_status))
                continue
            
            # Determine overall status for this sample
            all_r1_present = all(lane_status[lane]['R1']['exists'] for lane in existing_lanes)
            all_r2_present = all(lane_status[lane]['R2']['exists'] for lane in existing_lanes)
            
            if all_r1_present and all_r2_present:
                # Valid sample - both reads present in all existing lanes
                file_status[source]['valid'][sample_id] = {
                    'R1': {lane: lane_status[lane]['R1']['path'] for lane in existing_lanes},
                    'R2': {lane: lane_status[lane]['R2']['path'] for lane in existing_lanes},
                    'lane_status': lane_status
                }
            elif any(lane_status[lane]['R1']['exists'] for lane in existing_lanes) or \
                 any(lane_status[lane]['R2']['exists'] for lane in existing_lanes):
                # Incomplete - some files present but not all
                file_status[source]['incomplete'][sample_id] = {
                    'barcode': barcode,
                    'lane_status': lane_status
                }
            else:
                # Missing - no files found
                file_status[source]['missing'].append((sample_id, barcode, lane_status))
        
        # Report for this source
        n_valid = len(file_status[source]['valid'])
        n_incomplete = len(file_status[source]['incomplete'])
        n_missing = len(file_status[source]['missing'])
        
        status_parts = []
        if n_valid > 0:
            status_parts.append(f"[green]{n_valid} complete pairs[/green]")
        if n_incomplete > 0:
            status_parts.append(f"[yellow]{n_incomplete} incomplete[/yellow]")
        if n_missing > 0:
            status_parts.append(f"[red]{n_missing} missing[/red]")
        
        status_msg = f"  {source}: " + "  ".join(status_parts) if status_parts else f"  {source}: [dim]no samples[/dim]"
        console.print(status_msg)
    
    console.print()
    
    # ═══════════════════════════════════════════════════════════════
    # STEP 4: SUMMARY OF BROKEN SAMPLES
    # ═══════════════════════════════════════════════════════════════
    total_incomplete = sum(len(file_status[s]['incomplete']) for s in file_status)
    total_missing = sum(len(file_status[s]['missing']) for s in file_status)
    
    if total_incomplete > 0 or total_missing > 0:
        console.print(Panel(
            f"[yellow bold]⚠ ISSUES DETECTED[/yellow bold]\n\n"
            f"[yellow]Incomplete samples:[/yellow] {total_incomplete}\n"
            f"[red]Missing samples:[/red] {total_missing}",
            title="⚠ Sample Status Summary",
            border_style="yellow",
            box=box.DOUBLE
        ))
        
        # Create summary table of all broken samples
        summary_table = Table(
            box=box.ROUNDED,
            show_header=True,
            header_style="bold cyan",
            border_style="yellow",
            title="Samples with Issues"
        )
        summary_table.add_column("Sample ID", style="cyan", no_wrap=True)
        summary_table.add_column("Barcode", style="yellow", no_wrap=True, width=10)
        summary_table.add_column("Source", style="magenta", width=12)
        summary_table.add_column("Status", style="bold", width=12)
        summary_table.add_column("Issue Summary", style="dim")
        
        # Collect all broken samples
        for source in file_status:
            source_display = f"MERGED" if source == 'MERGED' else source
            
            # Add incomplete samples
            for sample_id, info in sorted(file_status[source]['incomplete'].items()):
                lane_status = info['lane_status']
                barcode = info['barcode']
                
                # Count missing files
                missing_count = 0
                for lane in lane_status:
                    if not lane_status[lane].get('dir_exists', False):
                        missing_count += 2
                    else:
                        if not lane_status[lane]['R1']['exists']:
                            missing_count += 1
                        if not lane_status[lane]['R2']['exists']:
                            missing_count += 1
                
                total_expected = len(lane_status) * 2
                found_count = total_expected - missing_count
                
                summary_table.add_row(
                    sample_id,
                    str(barcode),
                    source_display,
                    "[yellow]Incomplete[/yellow]",
                    f"{found_count}/{total_expected} files found"
                )
            
            # Add missing samples
            for sample_id, barcode, lane_status in sorted(file_status[source]['missing']):
                total_expected = len(lane_status) * 2
                
                summary_table.add_row(
                    sample_id,
                    str(barcode),
                    source_display,
                    "[red]Missing[/red]",
                    f"0/{total_expected} files found"
                )
        
        console.print(summary_table)
        console.print()
    
    # ═══════════════════════════════════════════════════════════════
    # STEP 5: DETAILED BREAKDOWN FOR BROKEN SAMPLES
    # ═══════════════════════════════════════════════════════════════
    
    # ─────────────────────────────────────────────────────────────
    # 5a. Incomplete samples (missing some reads/lanes)
    # ─────────────────────────────────────────────────────────────
    if total_incomplete > 0:
        console.print(Panel(
            "[cyan bold]Detailed Breakdown: Incomplete Samples[/cyan bold]\n\n"
            "[dim]Showing which specific files are present/missing for each sample[/dim]",
            border_style="cyan",
            box=box.ROUNDED
        ))
        
        for source in file_status:
            if file_status[source]['incomplete']:
                source_display = f"MERGED ({', '.join(merged_lanes)})" if source == 'MERGED' else source
                console.print(f"\n[bold yellow]Source: {source_display}[/bold yellow]")
                
                for sample_id, info in sorted(file_status[source]['incomplete'].items()):
                    lane_status = info['lane_status']
                    barcode = info['barcode']
                    
                    # Create a table showing status per lane
                    table = Table(
                        title=f"Sample: {sample_id} (Barcode: {barcode})",
                        box=box.ROUNDED,
                        show_header=True,
                        header_style="bold cyan",
                        border_style="yellow",
                        show_lines=True
                    )
                    table.add_column("Lane", style="cyan", width=6)
                    table.add_column("R1 Status", justify="center", width=10)
                    table.add_column("R1 Size", justify="right", width=10)
                    table.add_column("R2 Status", justify="center", width=10)
                    table.add_column("R2 Size", justify="right", width=10)
                    table.add_column("File Path", style="dim", overflow="fold")
                    
                    for lane in sorted(lane_status.keys()):
                        if not lane_status[lane].get('dir_exists', False):
                            lane_dir = os.path.join(base_dir, 'MERGED_LANES' if source == 'MERGED' else lane)
                            table.add_row(
                                lane,
                                "[red]✗[/red]",
                                "—",
                                "[red]✗[/red]",
                                "—",
                                f"[red]Directory not found: {lane_dir}[/red]"
                            )
                            continue
                        
                        r1_info = lane_status[lane]['R1']
                        r2_info = lane_status[lane]['R2']
                        
                        r1_status = "[green]✓[/green]" if r1_info['exists'] else "[red]✗[/red]"
                        r2_status = "[green]✓[/green]" if r2_info['exists'] else "[red]✗[/red]"
                        
                        # Add warning for small files
                        if r1_info.get('is_small'):
                            r1_status = "[yellow]⚠[/yellow]"
                        if r2_info.get('is_small'):
                            r2_status = "[yellow]⚠[/yellow]"
                        
                        r1_size = r1_info['size'] if r1_info['size'] else "—"
                        r2_size = r2_info['size'] if r2_info['size'] else "—"
                        
                        # Show path of missing file
                        missing_paths = []
                        if not r1_info['exists']:
                            missing_paths.append(f"R1: {r1_info['path']}")
                        if not r2_info['exists']:
                            missing_paths.append(f"R2: {r2_info['path']}")
                        
                        path_display = " | ".join(missing_paths) if missing_paths else os.path.dirname(r1_info['path'])
                        
                        table.add_row(lane, r1_status, r1_size, r2_status, r2_size, path_display)
                    
                    console.print(table)
                    console.print()
    
    # ─────────────────────────────────────────────────────────────
    # 5b. Completely missing samples (no files found at all)
    # ─────────────────────────────────────────────────────────────
    if total_missing > 0:
        if on_missing == 'error':
            console.print(Panel(
                "[red bold]❌ ERROR: MISSING FILES DETECTED[/red bold]\n\n"
                f"Found [red bold]{total_missing}[/red bold] sample(s) with NO files.\n"
                "[yellow]Execution stopped because on_missing='error'[/yellow]",
                title="❌ Stopping Execution",
                border_style="red",
                box=box.DOUBLE
            ))
        else:
            console.print(Panel(
                "[cyan bold]Detailed Breakdown: Completely Missing Samples[/cyan bold]\n\n"
                "[dim]Showing expected file paths for samples with no files found[/dim]",
                border_style="cyan",
                box=box.ROUNDED
            ))
        
        for source in file_status:
            if file_status[source]['missing']:
                source_display = f"MERGED ({', '.join(merged_lanes)})" if source == 'MERGED' else source
                console.print(f"\n[bold red]Source: {source_display}[/bold red]")
                
                for sample_id, barcode, lane_status in sorted(file_status[source]['missing']):
                    # Create a table showing what was expected in each lane
                    table = Table(
                        title=f"Sample: {sample_id} (Barcode: {barcode})",
                        box=box.ROUNDED,
                        show_header=True,
                        header_style="bold cyan",
                        border_style="red",
                        show_lines=True
                    )
                    table.add_column("Lane", style="cyan", width=6)
                    table.add_column("R1 Status", justify="center", width=10)
                    table.add_column("R2 Status", justify="center", width=10)
                    table.add_column("Expected File Paths", style="dim", overflow="fold")
                    
                    for lane in sorted(lane_status.keys()):
                        if not lane_status[lane].get('dir_exists', False):
                            lane_dir = os.path.join(base_dir, 'MERGED_LANES' if source == 'MERGED' else lane)
                            table.add_row(
                                lane,
                                "[red]✗[/red]",
                                "[red]✗[/red]",
                                f"[red]Directory not found: {lane_dir}[/red]"
                            )
                            continue
                        
                        r1_info = lane_status[lane]['R1']
                        r2_info = lane_status[lane]['R2']
                        
                        table.add_row(
                            lane,
                            "[red]✗[/red]",
                            "[red]✗[/red]",
                            f"R1: {r1_info['path']}\nR2: {r2_info['path']}"
                        )
                    
                    console.print(table)
                    console.print()
        
        if on_missing == 'error':
            return None
    
    # ═══════════════════════════════════════════════════════════════
    # STEP 6: CREATE SYMLINKS FOR VALID SAMPLES
    # ═══════════════════════════════════════════════════════════════
    if not dry_run:
        os.makedirs(target_dir, exist_ok=True)
        console.print(f"[green]✓[/green] Target directory ready: [cyan]{target_dir}[/cyan]\n")
    
    stats = {
        'total_samples': 0,
        'total_links': 0,
        'created': 0,
        'failed': 0,
        'skipped_incomplete': total_incomplete,
        'skipped_missing': total_missing,
        'small_files_warning': total_small_files,
        'by_source': {}
    }
    
    for source in file_status:
        valid_samples = file_status[source]['valid']
        
        if not valid_samples:
            continue
        
        # Header for this source
        source_display = f"MERGED ({', '.join(merged_lanes)})" if source == 'MERGED' else source
        
        if dry_run:
            console.print(Panel(
                f"[yellow bold]Symlinks that would be created for {source_display}[/yellow bold]\n"
                f"[dim]{len(valid_samples)} samples, {len(valid_samples) * 2} total links[/dim]",
                border_style="yellow",
                box=box.ROUNDED
            ))
        else:
            console.print(Panel(
                f"[green bold]Creating symlinks for {source_display}[/green bold]\n"
                f"[dim]{len(valid_samples)} samples, {len(valid_samples) * 2} total links[/dim]",
                border_style="green",
                box=box.ROUNDED
            ))
        
        table = Table(
            box=box.ROUNDED,
            show_header=True,
            header_style="bold cyan",
            show_lines=False,
            border_style="cyan" if not dry_run else "yellow"
        )
        table.add_column("#", justify="right", style="dim", width=4)
        table.add_column("Sample ID", style="cyan bold", no_wrap=True, width=25)
        table.add_column("R1", justify="center", width=6)
        table.add_column("R2", justify="center", width=6)
        table.add_column("Status", style="dim", width=20)
        
        source_stats = {'samples': len(valid_samples), 'created': 0, 'failed': 0}
        
        for idx, sample_id in enumerate(sorted(valid_samples.keys()), 1):
            info = valid_samples[sample_id]
            statuses = []
            status_msg = []
            
            for read in ['R1', 'R2']:
                # Get the source file path (use first lane from the dict)
                lanes_in_info = list(info[read].keys())
                source_file = info[read][lanes_in_info[0]]
                
                target_file = os.path.join(target_dir, f"{sample_id}_{read}.fastq.gz")
                
                if dry_run:
                    statuses.append("[yellow]○[/yellow]")
                    status_msg.append(f"{read}:pending")
                else:
                    try:
                        if os.path.islink(target_file) or os.path.exists(target_file):
                            os.remove(target_file)
                        os.symlink(source_file, target_file)
                        statuses.append("[green]✓[/green]")
                        status_msg.append(f"{read}:created")
                        source_stats['created'] += 1
                        stats['created'] += 1
                    except Exception as e:
                        statuses.append("[red]✗[/red]")
                        status_msg.append(f"{read}:failed")
                        source_stats['failed'] += 1
                        stats['failed'] += 1
                        console.print(f"[red]  Error: {sample_id} {read}: {e}[/red]")
            
            overall_status = " | ".join(status_msg)
            
            table.add_row(
                str(idx),
                sample_id,
                statuses[0],
                statuses[1],
                overall_status
            )
        
        console.print(table)
        
        # Source summary with nice box
        if dry_run:
            summary_box = Panel(
                f"[yellow]Would process {source_stats['samples']} samples ({source_stats['samples'] * 2} links)[/yellow]",
                border_style="yellow",
                box=box.ROUNDED,
                padding=(0, 2)
            )
        else:
            if source_stats['failed'] == 0:
                summary_box = Panel(
                    f"[green]✓ Successfully created {source_stats['created']} links for {source_stats['samples']} samples[/green]",
                    border_style="green",
                    box=box.ROUNDED,
                    padding=(0, 2)
                )
            else:
                summary_box = Panel(
                    f"[yellow]Created: {source_stats['created']} | Failed: {source_stats['failed']}[/yellow]",
                    border_style="yellow",
                    box=box.ROUNDED,
                    padding=(0, 2)
                )
        
        console.print(summary_box)
        console.print()
        
        stats['by_source'][source] = source_stats
        stats['total_samples'] += source_stats['samples']
        stats['total_links'] += source_stats['samples'] * 2
    
    # ═══════════════════════════════════════════════════════════════
    # STEP 7: FINAL SUMMARY
    # ═══════════════════════════════════════════════════════════════
    summary_table = Table(
        box=box.ROUNDED,
        show_header=False,
        border_style="cyan"
    )
    summary_table.add_column("", style="bold cyan", width=25)
    summary_table.add_column("", style="yellow", justify="right")
    
    if dry_run:
        summary_table.add_row("Mode", "[yellow]DRY RUN[/yellow]")
        summary_table.add_row("Would create links", str(stats['total_links']))
        summary_table.add_row("For samples", str(stats['total_samples']))
    else:
        summary_table.add_row("Mode", "[green]EXECUTION[/green]")
        summary_table.add_row("Links created", f"[green]{stats['created']}[/green]")
        summary_table.add_row("Links failed", f"[red]{stats['failed']}[/red]" if stats['failed'] > 0 else "0")
        summary_table.add_row("Samples processed", str(stats['total_samples']))
    
    if stats['skipped_incomplete'] > 0:
        summary_table.add_row("Skipped (incomplete)", f"[yellow]{stats['skipped_incomplete']}[/yellow]")
    
    if stats['skipped_missing'] > 0:
        summary_table.add_row("Skipped (missing)", f"[yellow]{stats['skipped_missing']}[/yellow]")
    
    if stats['small_files_warning'] > 0:
        summary_table.add_row("Small files (<" + str(min_file_size_mb) + " MB)", f"[yellow]{stats['small_files_warning']}[/yellow]")
    
    summary_table.add_row("Sources processed", str(len(stats['by_source'])))
    
    # Determine overall status
    if dry_run:
        border_color = "yellow"
        title = "📊 Preview Complete"
        footer_msg = "[dim]Set dry_run=False to execute[/dim]"
    elif stats['failed'] == 0 and stats['skipped_incomplete'] == 0 and stats['skipped_missing'] == 0 and stats['small_files_warning'] == 0:
        border_color = "green"
        title = "✅ Success"
        footer_msg = "[green]All symlinks created successfully![/green]"
    elif stats['failed'] > 0:
        border_color = "red"
        title = "⚠ Completed with Errors"
        footer_msg = f"[yellow]{stats['failed']} link(s) failed - check errors above[/yellow]"
    else:
        border_color = "yellow"
        title = "✅ Completed with Warnings"
        warnings = []
        if stats['skipped_incomplete'] > 0:
            warnings.append(f"{stats['skipped_incomplete']} incomplete")
        if stats['skipped_missing'] > 0:
            warnings.append(f"{stats['skipped_missing']} missing")
        if stats['small_files_warning'] > 0:
            warnings.append(f"{stats['small_files_warning']} small files")
        footer_msg = f"[yellow]Warnings: {', '.join(warnings)}[/yellow]"
    
    console.print(Panel(
        summary_table,
        title=title,
        border_style=border_color,
        box=box.DOUBLE
    ))
    
    console.print(f"\n{footer_msg}\n")
    
    return stats
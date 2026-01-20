import re
import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box


def rename_lane_columns(df, verbose=True):
    """
    Intelligently rename columns to include lane suffixes.
    
    Handles both naming patterns:
    Pattern 1: Barcode, BX, BX.1, BX.2
    Pattern 2: Barcode, Barcode.1, Barcode.2, Barcode.3
    
    Renamed to standardized format:
    - Barcode_L01, ID_L01, Target_L01, Progect_L01, %run_L01
    - Barcode_L02, ID_L02, Target_L02, Progect_L02, %run_L02
    - Barcode_L03, ID_L03, Target_L03, Progect_L03, %run_L03
    - Barcode_L04, ID_L04, Target_L04, Progect_L04, %run_L04
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with multi-lane sequencing data
    verbose : bool
        If True, print detailed information about detected patterns
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with renamed columns
    dict
        Mapping information for audit trail
    """
    console = Console()
    df = df.copy()
    
    if verbose:
        console.print(Panel(
            "[cyan]Auto-detecting lane column patterns...[/cyan]",
            title="Lane Column Renamer",
            box=box.ROUNDED
        ))
    
    # Base field names we expect (in standardized form)
    base_fields = ['Barcode', 'ID', 'Target', 'Progect', '%run']
    
    # Alternative names that might be used
    barcode_variants = ['Barcode', 'BX', 'barcode', 'bx']
    
    # Build rename mapping
    rename_map = {}
    detected_patterns = {}
    
    # For each base field, find all its variants in the dataframe
    for base_field in base_fields:
        # Get all possible names for this field
        if base_field == 'Barcode':
            search_names = barcode_variants
        else:
            search_names = [base_field]
        
        # Find columns matching this field
        field_columns = []

        # Special handling for Barcode field: check if "Barcode" column exists
        has_barcode_column = False
        if base_field == 'Barcode':
            has_barcode_column = 'Barcode' in [c.strip() for c in df.columns]

        for col in df.columns:
            col_stripped = col.strip()
            
            # Check if this column matches any of our search names
            for search_name in search_names:
                # Exact match handling
                if col_stripped == search_name:
                    # Special case: if we have "Barcode" column, then "BX" becomes L02
                    if base_field == 'Barcode' and search_name == 'BX' and has_barcode_column:
                        field_columns.append((col, 2, 'exact_shifted'))
                    else:
                        field_columns.append((col, 1, 'exact'))
                    break
                
                # Pattern: name.N
                pattern = rf'^{re.escape(search_name)}\.(\d+)\s*$'
                match = re.match(pattern, col_stripped)
                if match:
                    suffix_num = int(match.group(1))
                    # If BX and Barcode exists, shift BX.N lanes by 1
                    if base_field == 'Barcode' and search_name in ['BX', 'bx'] and has_barcode_column:
                        lane_num = suffix_num + 2  # BX.1 → L03, BX.2 → L04
                    else:
                        lane_num = suffix_num + 1  # Normal: BX.1 → L02
                    field_columns.append((col, lane_num, f'.{suffix_num}'))
                    break
        
        # Map to standardized names
        if field_columns:
            detected_patterns[base_field] = []
            for col, lane_num, suffix in field_columns:
                new_name = f'{base_field}_L{lane_num:02d}'
                rename_map[col] = new_name
                detected_patterns[base_field].append({
                    'original': col,
                    'lane': lane_num,
                    'suffix': suffix,
                    'new_name': new_name
                })
    
    # Show detected patterns
    if verbose and detected_patterns:
        for base_field, patterns in detected_patterns.items():
            table = Table(
                title=f"Field: {base_field}",
                box=box.SIMPLE,
                show_header=True,
                header_style="bold cyan"
            )
            table.add_column("Original Column", style="yellow")
            table.add_column("Lane", justify="center", style="cyan")
            table.add_column("→", justify="center", style="dim")
            table.add_column("New Column", style="green")
            
            for pattern in sorted(patterns, key=lambda x: x['lane']):
                table.add_row(
                    f"'{pattern['original']}'",
                    f"L{pattern['lane']:02d}",
                    "→",
                    pattern['new_name']
                )
            
            console.print(table)
    
    # Apply renaming
    df_renamed = df.rename(columns=rename_map)
    
    # Summary
    if verbose:
        lanes_found = set()
        for patterns in detected_patterns.values():
            lanes_found.update(p['lane'] for p in patterns)
        
        lanes_str = ', '.join(f"L{l:02d}" for l in sorted(lanes_found))
        
        summary_table = Table(box=box.ROUNDED, show_header=False, border_style="green")
        summary_table.add_column("", style="bold")
        summary_table.add_column("", style="cyan")
        
        summary_table.add_row("Fields detected:", ', '.join(sorted(detected_patterns.keys())))
        summary_table.add_row("Lanes found:", lanes_str)
        summary_table.add_row("Columns renamed:", str(len(rename_map)))
        
        console.print(Panel(
            summary_table,
            title="✓ Renaming Complete",
            border_style="green",
            box=box.ROUNDED
        ))
    
    # Build audit trail
    audit = {
        'renamed_columns': rename_map,
        'detected_patterns': detected_patterns,
        'lanes_found': sorted(set(p['lane'] for patterns in detected_patterns.values() for p in patterns))
    }
    
    return df_renamed, audit


def extract_lane_columns(df_renamed, lanes=['L01', 'L02', 'L03', 'L04'], verbose=True):
    """
    Extract only the columns for specified lanes.
    
    Parameters:
    -----------
    df_renamed : pd.DataFrame
        DataFrame with renamed columns (output from rename_lane_columns)
    lanes : list
        List of lanes to extract (e.g., ['L02', 'L03'])
    verbose : bool
        If True, print information about extracted columns
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with only the specified lane columns
    """
    console = Console()
    
    columns_to_keep = []
    found_by_lane = {}
    
    for lane in lanes:
        lane_cols = [col for col in df_renamed.columns if col.endswith(f'_{lane}')]
        columns_to_keep.extend(lane_cols)
        found_by_lane[lane] = lane_cols
    
    if verbose:
        table = Table(
            title="Extracting Lane Columns",
            box=box.ROUNDED,
            show_header=True,
            header_style="bold cyan"
        )
        table.add_column("Lane", style="cyan")
        table.add_column("Columns Found", justify="right", style="yellow")
        table.add_column("Column Names", style="dim")
        
        for lane in lanes:
            cols = found_by_lane.get(lane, [])
            col_names = ', '.join(cols) if cols else "[red]none[/red]"
            table.add_row(lane, str(len(cols)), col_names)
        
        console.print(table)
        
        if columns_to_keep:
            console.print(f"\n[green]✓ Extracting {len(columns_to_keep)} columns from {len(lanes)} lane(s)[/green]")
        else:
            console.print(f"\n[yellow]⚠ No columns found matching the specified lanes[/yellow]")
    
    return df_renamed[columns_to_keep]
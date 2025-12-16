

import re
import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box



def validate_and_clean_ids(df, auto_fix_cyrillic=True, auto_fix_spaces=True, auto_strip=True):
    """
    Validate ID columns for invalid characters and optionally fix common issues.
    
    Valid characters: Latin letters (a-z, A-Z), numbers (0-9), underscore (_), hyphen (-)
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with ID columns
    auto_fix_cyrillic : bool
        Automatically replace Cyrillic lookalikes with Latin equivalents
    auto_fix_spaces : bool
        Automatically replace spaces with underscores
    auto_strip : bool
        Automatically strip leading/trailing whitespace
    
    Returns:
    --------
    pd.DataFrame
        Cleaned DataFrame
    dict
        Report of issues found and fixed
    """
    console = Console()
    df = df.copy()
    
    # Cyrillic to Latin mapping
    cyrillic_to_latin = {
        'А': 'A', 'В': 'B', 'С': 'C', 'Е': 'E', 'Н': 'H',
        'К': 'K', 'М': 'M', 'О': 'O', 'Р': 'P', 'Т': 'T',
        'Х': 'X', 'У': 'Y',
        'а': 'a', 'в': 'b', 'е': 'e', 'о': 'o', 'р': 'p', 
        'с': 'c', 'у': 'y', 'х': 'x', 'к': 'k', 'н': 'h'
    }
    
    # Find all ID columns
    id_columns = [col for col in df.columns if col.startswith('ID')]
    
    if not id_columns:
        console.print("[yellow]No ID columns found in DataFrame[/yellow]")
        return df, {}
    
    console.print(Panel(
        f"[cyan]Validating {len(id_columns)} ID column(s): {', '.join(id_columns)}[/cyan]",
        box=box.ROUNDED
    ))
    
    # Pattern for valid characters: Latin letters, numbers, underscore, hyphen
    valid_pattern = re.compile(r'^[a-zA-Z0-9_-]+$')
    
    report = {
        'cyrillic_fixed': {},
        'spaces_fixed': {},
        'whitespace_stripped': {},
        'invalid_chars': {},
        'total_issues': 0
    }
    
    for col in id_columns:
        if col not in df.columns:
            continue
        
        console.print(f"\n[bold cyan]Checking column: {col}[/bold cyan]")
        
        # Convert to string and filter out NaN values
        col_values = df[col].astype(str)
        non_null_mask = df[col].notna()
        
        issues_in_column = []
        
        # Check each value
        for idx, value in col_values[non_null_mask].items():
            if value == 'nan':
                continue
            
            original_value = value
            modified = False
            
            # 0. Check and strip leading/trailing whitespace
            if value != value.strip():
                if col not in report['whitespace_stripped']:
                    report['whitespace_stripped'][col] = []
                report['whitespace_stripped'][col].append({
                    'row': idx,
                    'original': value,
                    'fixed': value.strip() if auto_strip else None
                })
                if auto_strip:
                    value = value.strip()
                    modified = True
                report['total_issues'] += 1
            
            # 1. Check and fix Cyrillic characters
            cyrillic_found = []
            for cyrillic_char in cyrillic_to_latin.keys():
                if cyrillic_char in value:
                    cyrillic_found.append(cyrillic_char)
                    if auto_fix_cyrillic:
                        value = value.replace(cyrillic_char, cyrillic_to_latin[cyrillic_char])
                        modified = True
            
            if cyrillic_found:
                if col not in report['cyrillic_fixed']:
                    report['cyrillic_fixed'][col] = []
                report['cyrillic_fixed'][col].append({
                    'row': idx,
                    'original': original_value,
                    'fixed': value if auto_fix_cyrillic else None,
                    'chars': cyrillic_found
                })
                report['total_issues'] += 1
            
            # 2. Check and fix spaces
            if ' ' in value:
                if col not in report['spaces_fixed']:
                    report['spaces_fixed'][col] = []
                report['spaces_fixed'][col].append({
                    'row': idx,
                    'original': original_value if not modified else value,
                    'fixed': value.replace(' ', '_') if auto_fix_spaces else None
                })
                if auto_fix_spaces:
                    value = value.replace(' ', '_')
                    modified = True
                report['total_issues'] += 1
            
            # 3. Check for other invalid characters
            if not valid_pattern.match(value):
                invalid_chars = set()
                for char in value:
                    if not re.match(r'[a-zA-Z0-9_-]', char):
                        invalid_chars.add(char)
                
                if invalid_chars:
                    if col not in report['invalid_chars']:
                        report['invalid_chars'][col] = []
                    report['invalid_chars'][col].append({
                        'row': idx,
                        'value': value,
                        'invalid_chars': list(invalid_chars)
                    })
                    issues_in_column.append((idx, value, invalid_chars))
                    report['total_issues'] += 1
            
            # Update the dataframe if modified
            if modified:
                df.at[idx, col] = value
        
        # Display results for this column
        if col in report['cyrillic_fixed']:
            table = Table(title=f"Cyrillic Characters {'Fixed' if auto_fix_cyrillic else 'Found'}", 
                         box=box.ROUNDED, show_lines=True)
            table.add_column("Row", style="cyan", width=8)
            table.add_column("Original", style="yellow")
            table.add_column("Cyrillic Chars", style="red")
            if auto_fix_cyrillic:
                table.add_column("Fixed", style="green")
            
            for item in report['cyrillic_fixed'][col][:10]:  # Show first 10
                row_str = str(item['row'])
                chars_str = ', '.join(item['chars'])
                if auto_fix_cyrillic:
                    table.add_row(row_str, item['original'], chars_str, item['fixed'])
                else:
                    table.add_row(row_str, item['original'], chars_str)
            
            if len(report['cyrillic_fixed'][col]) > 10:
                console.print(table)
                console.print(f"[dim]... and {len(report['cyrillic_fixed'][col]) - 10} more[/dim]")
            else:
                console.print(table)
        
        if col in report['whitespace_stripped']:
            table = Table(title=f"Leading/Trailing Whitespace {'Fixed' if auto_strip else 'Found'}", 
                         box=box.ROUNDED, show_lines=True)
            table.add_column("Row", style="cyan", width=8)
            table.add_column("Original", style="yellow")
            if auto_strip:
                table.add_column("Fixed", style="green")
            
            for item in report['whitespace_stripped'][col][:10]:
                # Show whitespace visually with special characters
                original_display = item['original'].replace(' ', '·')
                if auto_strip:
                    fixed_display = item['fixed'].replace(' ', '·') if item['fixed'] else ''
                    table.add_row(str(item['row']), f"'{original_display}'", f"'{fixed_display}'")
                else:
                    table.add_row(str(item['row']), f"'{original_display}'")
            
            if len(report['whitespace_stripped'][col]) > 10:
                console.print(table)
                console.print(f"[dim]... and {len(report['whitespace_stripped'][col]) - 10} more[/dim]")
            else:
                console.print(table)
        
        if col in report['spaces_fixed']:
            table = Table(title=f"Spaces {'Fixed' if auto_fix_spaces else 'Found'}", 
                         box=box.ROUNDED, show_lines=True)
            table.add_column("Row", style="cyan", width=8)
            table.add_column("Original", style="yellow")
            if auto_fix_spaces:
                table.add_column("Fixed", style="green")
            
            for item in report['spaces_fixed'][col][:10]:
                if auto_fix_spaces:
                    table.add_row(str(item['row']), item['original'], item['fixed'])
                else:
                    table.add_row(str(item['row']), item['original'])
            
            if len(report['spaces_fixed'][col]) > 10:
                console.print(table)
                console.print(f"[dim]... and {len(report['spaces_fixed'][col]) - 10} more[/dim]")
            else:
                console.print(table)
        
        if col in report['invalid_chars']:
            table = Table(title="⚠️  Invalid Characters Detected", 
                         box=box.ROUNDED, show_lines=True, 
                         border_style="red")
            table.add_column("Row", style="cyan", width=8)
            table.add_column("Value", style="yellow")
            table.add_column("Invalid Characters", style="red bold")
            
            for item in report['invalid_chars'][col][:10]:
                chars_display = ', '.join([f"'{c}' (U+{ord(c):04X})" for c in item['invalid_chars']])
                table.add_row(str(item['row']), item['value'], chars_display)
            
            if len(report['invalid_chars'][col]) > 10:
                console.print(table)
                console.print(f"[dim]... and {len(report['invalid_chars'][col]) - 10} more[/dim]")
            else:
                console.print(table)
        
        # Column summary
        if (col not in report['cyrillic_fixed'] and 
            col not in report['whitespace_stripped'] and
            col not in report['spaces_fixed'] and 
            col not in report['invalid_chars']):
            console.print(f"[green]✓ Column '{col}' is clean (all IDs use valid characters)[/green]")
    
    # Overall summary
    console.print("\n")
    summary_table = Table(title="Validation Summary", box=box.DOUBLE, border_style="cyan")
    summary_table.add_column("Issue Type", style="bold")
    summary_table.add_column("Count", justify="right", style="yellow")
    summary_table.add_column("Status", style="green")
    
    cyrillic_count = sum(len(v) for v in report['cyrillic_fixed'].values())
    whitespace_count = sum(len(v) for v in report['whitespace_stripped'].values())
    spaces_count = sum(len(v) for v in report['spaces_fixed'].values())
    invalid_count = sum(len(v) for v in report['invalid_chars'].values())
    
    if cyrillic_count > 0:
        status = "✓ Auto-fixed" if auto_fix_cyrillic else "⚠ Needs fixing"
        summary_table.add_row("Cyrillic characters", str(cyrillic_count), status)
    
    if whitespace_count > 0:
        status = "✓ Auto-fixed" if auto_strip else "⚠ Needs fixing"
        summary_table.add_row("Leading/trailing spaces", str(whitespace_count), status)
    
    if spaces_count > 0:
        status = "✓ Auto-fixed" if auto_fix_spaces else "⚠ Needs fixing"
        summary_table.add_row("Spaces in IDs", str(spaces_count), status)
    
    if invalid_count > 0:
        summary_table.add_row("Other invalid chars", str(invalid_count), "⚠ Manual fix required")
    
    if report['total_issues'] == 0:
        summary_table.add_row("All IDs", "—", "✓ Valid")
    
    console.print(summary_table)
    
    # Warning message if there are unfixed issues
    if invalid_count > 0:
        console.print(Panel(
            "[bold red]⚠️  WARNING[/bold red]\n\n"
            "Some IDs contain invalid characters that cannot be auto-fixed.\n"
            "Valid characters: Latin letters (a-z, A-Z), numbers (0-9), underscore (_), hyphen (-)\n\n"
            "[yellow]Please manually correct these IDs in your source data.[/yellow]",
            border_style="red",
            box=box.DOUBLE
        ))
    
    return df, report


def normalize_barcode_ids(df, auto_normalize=True, verbose=True):
    """
    Normalize barcode IDs by extracting numeric values and removing leading zeros.
    
    Converts:
    - SP001 → 1
    - SP020 → 20
    - BC05 → 5
    - 001 → 1
    - 66 → 66 (already normalized)
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with Barcode columns
    auto_normalize : bool
        If True, automatically normalize all barcodes
        If False, only report what would be changed
    verbose : bool
        If True, print detailed information about changes
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with normalized barcodes
    dict
        Report of changes made
    """
    console = Console()
    df = df.copy()
    
    # Find all Barcode columns
    barcode_columns = [col for col in df.columns if col.startswith('Barcode')]
    
    if not barcode_columns:
        if verbose:
            console.print("[yellow]No Barcode columns found in DataFrame[/yellow]")
        return df, {}
    
    if verbose:
        console.print(Panel(
            f"[cyan]Normalizing {len(barcode_columns)} Barcode column(s): {', '.join(barcode_columns)}[/cyan]",
            title="Barcode Normalizer",
            box=box.ROUNDED
        ))
    
    # Pattern to extract numbers from barcodes
    # Matches optional prefix followed by digits
    number_pattern = re.compile(r'^[A-Za-z]*0*(\d+)$')
    
    report = {
        'normalized': {},  # column: [(original, normalized, row)]
        'already_clean': {},  # column: count
        'invalid': {},  # column: [(value, row)]
        'total_changes': 0
    }
    
    for col in barcode_columns:
        if col not in df.columns:
            continue
        
        if verbose:
            console.print(f"\n[bold cyan]Processing column: {col}[/bold cyan]")
        
        # Convert to string and filter out NaN
        col_values = df[col].astype(str)
        non_null_mask = df[col].notna()
        
        normalized_count = 0
        clean_count = 0
        invalid_count = 0
        
        changes = []
        invalid_values = []
        
        for idx, value in col_values[non_null_mask].items():
            if value == 'nan':
                continue
            
            original_value = value.strip()
            
            # Try to extract numeric part
            match = number_pattern.match(original_value)
            
            if match:
                # Extract the number without leading zeros
                normalized_value = match.group(1)
                
                # Check if it changed
                if original_value != normalized_value:
                    changes.append((original_value, normalized_value, idx))
                    normalized_count += 1
                    
                    if auto_normalize:
                        df.at[idx, col] = normalized_value
                else:
                    clean_count += 1
            else:
                # Couldn't parse as expected format
                invalid_values.append((original_value, idx))
                invalid_count += 1
        
        # Store results
        if changes:
            report['normalized'][col] = changes
            report['total_changes'] += len(changes)
        
        if clean_count > 0:
            report['already_clean'][col] = clean_count
        
        if invalid_values:
            report['invalid'][col] = invalid_values
        
        # Display results for this column
        if verbose:
            if changes:
                action = "Normalized" if auto_normalize else "Would normalize"
                table = Table(
                    title=f"{action} ({len(changes)} values)",
                    box=box.ROUNDED,
                    show_lines=True
                )
                table.add_column("Row", style="cyan", width=8)
                table.add_column("Original", style="yellow")
                table.add_column("→", justify="center", style="dim", width=3)
                table.add_column("Normalized", style="green")
                
                # Show first 10 examples
                for original, normalized, row_idx in changes[:10]:
                    table.add_row(str(row_idx), original, "→", normalized)
                
                console.print(table)
                
                if len(changes) > 10:
                    console.print(f"[dim]... and {len(changes) - 10} more[/dim]")
            
            if invalid_values:
                table = Table(
                    title="⚠️  Invalid Barcode Format",
                    box=box.ROUNDED,
                    show_lines=True,
                    border_style="red"
                )
                table.add_column("Row", style="cyan", width=8)
                table.add_column("Value", style="red bold")
                table.add_column("Issue", style="yellow")
                
                for value, row_idx in invalid_values[:10]:
                    table.add_row(
                        str(row_idx),
                        value,
                        "No numeric part found"
                    )
                
                console.print(table)
                
                if len(invalid_values) > 10:
                    console.print(f"[dim]... and {len(invalid_values) - 10} more[/dim]")
            
            # Column summary
            status_parts = []
            if clean_count > 0:
                status_parts.append(f"[green]{clean_count} already clean[/green]")
            if normalized_count > 0:
                status = "normalized" if auto_normalize else "to normalize"
                status_parts.append(f"[yellow]{normalized_count} {status}[/yellow]")
            if invalid_count > 0:
                status_parts.append(f"[red]{invalid_count} invalid[/red]")
            
            if status_parts:
                console.print("  " + "  ".join(status_parts))
            else:
                console.print("[green]✓ All values are clean[/green]")
    
    # Overall summary
    if verbose:
        console.print("\n")
        summary_table = Table(
            title="Normalization Summary",
            box=box.DOUBLE,
            border_style="cyan"
        )
        summary_table.add_column("Metric", style="bold")
        summary_table.add_column("Count", justify="right", style="yellow")
        summary_table.add_column("Status", style="green")
        
        total_normalized = sum(len(v) for v in report['normalized'].values())
        total_clean = sum(v for v in report['already_clean'].values())
        total_invalid = sum(len(v) for v in report['invalid'].values())
        
        if total_clean > 0:
            summary_table.add_row(
                "Already normalized",
                str(total_clean),
                "✓ No change needed"
            )
        
        if total_normalized > 0:
            status = "✓ Normalized" if auto_normalize else "⚠ Needs normalization"
            summary_table.add_row(
                "Values changed",
                str(total_normalized),
                status
            )
        
        if total_invalid > 0:
            summary_table.add_row(
                "Invalid format",
                str(total_invalid),
                "⚠ Manual review needed"
            )
        
        if total_normalized == 0 and total_invalid == 0:
            summary_table.add_row(
                "All barcodes",
                str(total_clean),
                "✓ Already normalized"
            )
        
        console.print(summary_table)
        
        # Warning for invalid barcodes
        if total_invalid > 0:
            console.print(Panel(
                "[bold red]⚠️  WARNING[/bold red]\n\n"
                f"Found {total_invalid} barcode(s) with invalid format.\n"
                "Expected format: Optional prefix (letters) + digits (e.g., BC001, SP20, 5)\n\n"
                "[yellow]Please review and correct these values manually.[/yellow]",
                border_style="red",
                box=box.DOUBLE
            ))
        
        # Info about what was done
        if not auto_normalize and total_normalized > 0:
            console.print(Panel(
                "[yellow]This was a preview.[/yellow]\n"
                f"Set auto_normalize=True to apply changes to {total_normalized} barcode(s).",
                border_style="yellow",
                box=box.ROUNDED
            ))
    
    return df, report
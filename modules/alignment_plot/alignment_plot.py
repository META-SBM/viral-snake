import pandas as pd
import glob
import re
from pathlib import Path
import click
import webbrowser
import tempfile
import os
import pickle
from typing import List, Dict, Optional
import plotly.graph_objects as go
import pysam

class Minimap2DataBrowser:
    """
    Browser for minimap2 alignment results with interactive exploration
    """
    
    def __init__(self, base_pattern: str = None):
        self.base_pattern = base_pattern
        self.df = None
        if base_pattern:
            self.load_data()
    
    def load_data(self) -> None:
        """Load and parse minimap2 log files"""
        click.echo("📁 Loading minimap2 data...")
        self.df = self.parse_minimap2_logs(self.base_pattern)
        click.echo(f"✅ Loaded {len(self.df)} alignment records")
    
    def parse_minimap2_logs(self, base_pattern: str) -> pd.DataFrame:
        """
        Parse all minimap2 log files and extract alignment statistics.
        """
        log_files = glob.glob(base_pattern, recursive=True)
        
        data = []
        
        for log_path in log_files:
            path_parts = Path(log_path).parts
            
            try:
                prefix_idx = [i for i, p in enumerate(path_parts) if 'contigs_formatted_minlen_' in p][0] - 1
                prefix = path_parts[prefix_idx]
                
                contigs_dir = [p for p in path_parts if 'contigs_formatted_minlen_' in p][0]
                min_len = contigs_dir.split('_')[-1]
                
                minimap_idx = [i for i, p in enumerate(path_parts) if p == 'minimap2'][0]
                reference_org = path_parts[minimap_idx + 1]
            except:
                click.echo(f"⚠️  Could not parse path: {log_path}")
                continue
            
            total_alignments = None
            mapped_contigs = None
            
            try:
                with open(log_path, 'r') as f:
                    lines = f.readlines()
                    
                for line in reversed(lines):
                    if 'Total alignments:' in line:
                        total_alignments = int(re.search(r'Total alignments:\s*(\d+)', line).group(1))
                    if 'Mapped contigs:' in line:
                        mapped_contigs = int(re.search(r'Mapped contigs:\s*(\d+)', line).group(1))
                    
                    if total_alignments is not None and mapped_contigs is not None:
                        break
            except Exception as e:
                click.echo(f"⚠️  Error reading {log_path}: {e}")
                continue
            
            png_path = str(Path(log_path).parent / 'bamsnap' / 'coverage.png')
            png_exists = Path(png_path).exists()
            
            data.append({
                'prefix': prefix,
                'min_len': min_len,
                'reference_org': reference_org,
                'total_alignments': total_alignments or 0,
                'mapped_contigs': mapped_contigs or 0,
                'log_path': log_path,
                'png_path': png_path,
                'png_exists': png_exists
            })
        
        df = pd.DataFrame(data)
        df = df.sort_values('mapped_contigs', ascending=False).reset_index(drop=True)
        
        return df
    
    def get_taxon_stats(self) -> Dict:
        """Get statistics by taxon"""
        if self.df is None or self.df.empty:
            return {}
        
        stats = self.df[self.df['mapped_contigs'] > 0].groupby('reference_org').agg({
            'mapped_contigs': ['count', 'sum', 'mean'],
            'total_alignments': 'sum'
        }).round(2)
        
        stats.columns = ['record_count', 'total_mapped_contigs', 'avg_mapped_contigs', 'total_alignments']
        stats = stats.sort_values('total_mapped_contigs', ascending=False)
        return stats.reset_index().to_dict('records')
    
    def search_taxon(self, taxon_query: str) -> pd.DataFrame:
        """Search for taxon by name (case insensitive partial match)"""
        if self.df is None:
            return pd.DataFrame()
        
        mask = self.df['reference_org'].str.contains(taxon_query, case=False, na=False)
        return self.df[mask].reset_index(drop=False)
    
    def get_record_by_index(self, index: int) -> Optional[Dict]:
        """Get record by index"""
        if self.df is None or index < 0 or index >= len(self.df):
            return None
        return self.df.iloc[index].to_dict()
    
    def plot_alignments_plotly(self, bam_path: str, reference: str = None, 
                             min_length: int = 0, max_contigs: int = None,
                             color_by: str = 'strand', title: str = None) -> Optional[go.Figure]:
        """
        Create interactive Plotly visualization of contig alignments.
        """
        bam_path = str(bam_path)
        
        if not Path(bam_path).exists():
            click.echo(f"❌ BAM file not found: {bam_path}")
            return None
        
        try:
            bam = pysam.AlignmentFile(bam_path, 'rb')
            
            refs_to_plot = [reference] if reference else list(bam.references)
            
            all_data = []
            ref_info = {}
            
            for ref_name in refs_to_plot:
                if ref_name not in bam.references:
                    continue
                    
                ref_len = bam.get_reference_length(ref_name)
                alignments = []
                
                for read in bam.fetch(ref_name):
                    if read.is_unmapped:
                        continue
                    
                    aln_len = read.reference_end - read.reference_start
                    if aln_len < min_length:
                        continue
                    
                    nm = read.get_tag('NM') if read.has_tag('NM') else None
                    if nm is not None and read.query_alignment_length > 0:
                        identity = (1 - nm / read.query_alignment_length) * 100
                    else:
                        matches = sum(length for op, length in read.cigartuples if op == 0)
                        identity = (matches / read.query_alignment_length * 100) if read.query_alignment_length > 0 else 0
                    
                    alignments.append({
                        'ref_name': ref_name,
                        'contig': read.query_name,
                        'start': read.reference_start,
                        'end': read.reference_end,
                        'length': aln_len,
                        'strand': '-' if read.is_reverse else '+',
                        'mapq': read.mapping_quality,
                        'identity': identity,
                        'query_start': read.query_alignment_start,
                        'query_end': read.query_alignment_end,
                        'query_length': read.query_length
                    })
                
                if max_contigs and len(alignments) > max_contigs:
                    alignments = sorted(alignments, key=lambda x: x['length'], reverse=True)[:max_contigs]
                
                ref_info[ref_name] = {
                    'length': ref_len, 
                    'n_alignments': len(alignments)
                }
                all_data.extend(alignments)
            
            if not all_data:
                click.echo("⚠️  No alignments found")
                return None
            
            df = pd.DataFrame(all_data)
            
            click.echo(f"📊 Total alignments: {len(df)}")
            click.echo(f"📊 Unique contigs: {df['contig'].nunique()}")
            click.echo(f"📊 References: {len(ref_info)}")
            
            fig = go.Figure()
            
            y_offset = 0
            annotations = []
            color_schemes = {'strand': {'+': '#3498db', '-': '#e74c3c'}}
            
            for ref_idx, (ref_name, info) in enumerate(ref_info.items()):
                ref_len = info['length']
                ref_data = df[df['ref_name'] == ref_name].copy()
                
                if len(ref_data) == 0:
                    continue
                
                ref_data = ref_data.sort_values('start').reset_index(drop=True)
                ref_data['y_pos'] = y_offset + 1
                
                for idx in range(len(ref_data)):
                    if idx > 0:
                        max_y = y_offset + 1
                        for prev_idx in range(idx):
                            prev_row = ref_data.iloc[prev_idx]
                            curr_row = ref_data.iloc[idx]
                            if not (curr_row['start'] > prev_row['end'] or curr_row['end'] < prev_row['start']):
                                max_y = max(max_y, prev_row['y_pos'] + 1)
                        ref_data.loc[ref_data.index[idx], 'y_pos'] = max_y
                
                max_y_pos = ref_data['y_pos'].max()
                
                fig.add_trace(go.Scatter(
                    x=[0, ref_len],
                    y=[y_offset, y_offset],
                    mode='lines',
                    line=dict(color='#2c3e50', width=8),
                    text=f"Reference: {ref_name}<br>Length: {ref_len:,} bp",
                    hoverinfo='text',
                    showlegend=False
                ))
                
                annotations.append(dict(
                    x=-ref_len * 0.02,
                    y=y_offset + (max_y_pos - y_offset) / 2,
                    text=f"<b>{ref_name}</b><br>{ref_len:,} bp",
                    showarrow=False,
                    xanchor='right',
                    font=dict(size=11)
                ))
                
                if color_by == 'strand':
                    for strand in ['+', '-']:
                        strand_data = ref_data[ref_data['strand'] == strand]
                        if len(strand_data) == 0:
                            continue
                        
                        color = color_schemes['strand'][strand]
                        
                        for _, row in strand_data.iterrows():
                            y_bottom = row['y_pos'] - 0.4
                            y_top = row['y_pos'] + 0.4
                            
                            hover_text = (
                                f"<b style='font-size:14px'>{row['contig']}</b><br>"
                                f"<br>"
                                f"<b>Reference Position:</b><br>"
                                f"  {row['start']:,} - {row['end']:,} bp<br>"
                                f"  Length: {row['length']:,} bp<br>"
                                f"<br>"
                                f"<b>Contig Info:</b><br>"
                                f"  Query region: {row['query_start']}-{row['query_end']} / {row['query_length']} bp<br>"
                                f"  Strand: <b>{row['strand']}</b><br>"
                                f"<br>"
                                f"<b>Alignment Quality:</b><br>"
                                f"  Identity: {row['identity']:.2f}%<br>"
                                f"  MAPQ: {row['mapq']}"
                            )
                            
                            fig.add_trace(go.Scatter(
                                x=[row['start'], row['end'], row['end'], row['start'], row['start']],
                                y=[y_bottom, y_bottom, y_top, y_top, y_bottom],
                                fill='toself',
                                fillcolor=color,
                                line=dict(color='#2c3e50', width=0.5),
                                mode='lines',
                                text=hover_text,
                                hoverinfo='text',
                                showlegend=False,
                                opacity=0.8
                            ))
                
                y_offset = max_y_pos + 2.5
            
            if color_by == 'strand':
                for strand, color in color_schemes['strand'].items():
                    fig.add_trace(go.Scatter(
                        x=[None], y=[None],
                        mode='markers',
                        marker=dict(size=10, color=color),
                        showlegend=True,
                        name=f'Strand {strand}'
                    ))
            
            total_height = max(400, int(y_offset * 70))
            
            fig.update_layout(
                title=title or "Contig Alignments",
                xaxis=dict(
                    title='Reference Position (bp)',
                    showgrid=True,
                    gridcolor='#ecf0f1',
                    tickformat=','
                ),
                yaxis=dict(
                    showticklabels=False,
                    showgrid=False
                ),
                annotations=annotations,
                height=total_height,
                hovermode='closest',
                plot_bgcolor='white',
                showlegend=(color_by == 'strand'),
                margin=dict(l=200, r=100, t=80, b=80)
            )
            
            bam.close()
            
            # Print alignment statistics
            click.echo(f"\n{'='*60}")
            click.echo("ALIGNMENT STATISTICS")
            click.echo(f"{'='*60}")
            for ref_name, info in ref_info.items():
                ref_data = df[df['ref_name'] == ref_name]
                click.echo(f"\n{ref_name}:")
                click.echo(f"  Alignments: {len(ref_data)}")
                click.echo(f"  Unique contigs: {ref_data['contig'].nunique()}")
                click.echo(f"  Avg identity: {ref_data['identity'].mean():.2f}%")
            
            click.echo(f"\n✅ Plot complete!")
            return fig
            
        except Exception as e:
            click.echo(f"❌ Error: {e}")
            import traceback
            traceback.print_exc()
            return None

    def save_to_file(self, filepath: str) -> None:
        """Save browser state to file"""
        with open(filepath, 'wb') as f:
            pickle.dump({
                'df': self.df,
                'base_pattern': self.base_pattern
            }, f)
        click.echo(f"✅ Browser state saved to: {filepath}")

    @classmethod
    def load_from_file(cls, filepath: str) -> 'Minimap2DataBrowser':
        """Load browser state from file"""
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        
        browser = cls()
        browser.df = data['df']
        browser.base_pattern = data['base_pattern']
        click.echo(f"✅ Browser state loaded from: {filepath}")
        return browser

# Global browser instance
_browser = None

def get_browser():
    """Get or create browser instance"""
    global _browser
    if _browser is None:
        # Try to load from default location
        default_file = "minimap2_browser_state.pkl"
        if os.path.exists(default_file):
            _browser = Minimap2DataBrowser.load_from_file(default_file)
        else:
            _browser = Minimap2DataBrowser()
    return _browser

def save_browser():
    """Save browser state to default location"""
    global _browser
    if _browser and _browser.df is not None:
        _browser.save_to_file("minimap2_browser_state.pkl")

# CLI Commands
@click.group()
def cli():
    """Minimap2 Alignment Results Browser"""
    pass

@cli.command()
@click.argument('base_pattern')
def load(base_pattern):
    """Load minimap2 data from pattern"""
    global _browser
    _browser = Minimap2DataBrowser(base_pattern)
    save_browser()
    click.echo(f"✅ Data loaded successfully with {len(_browser.df)} records")

@cli.command()
def stats():
    """Show taxon statistics"""
    browser = get_browser()
    if browser.df is None:
        click.echo("❌ No data loaded. Use 'load' command first.")
        return
    
    stats_data = browser.get_taxon_stats()
    if not stats_data:
        click.echo("❌ No statistics available")
        return
    
    click.echo("\n📊 TAXON STATISTICS")
    click.echo("=" * 80)
    click.echo(f"{'Taxon':<40} {'Records':<8} {'Mapped Contigs':<15} {'Avg Mapped':<12} {'Total Alignments':<15}")
    click.echo("-" * 80)
    
    for stat in stats_data:
        click.echo(f"{stat['reference_org']:<40} {stat['record_count']:<8} {stat['total_mapped_contigs']:<15} {stat['avg_mapped_contigs']:<12} {stat['total_alignments']:<15}")

@cli.command()
@click.argument('taxon_query')
def search(taxon_query):
    """Search for taxon by name"""
    browser = get_browser()
    if browser.df is None:
        click.echo("❌ No data loaded. Use 'load' command first.")
        return
    
    results = browser.search_taxon(taxon_query)
    if results.empty:
        click.echo(f"❌ No results found for '{taxon_query}'")
        return
    
    click.echo(f"\n🔍 SEARCH RESULTS FOR '{taxon_query.upper()}'")
    click.echo("=" * 100)
    click.echo(f"{'Index':<6} {'Taxon':<40} {'Min Len':<8} {'Mapped':<8} {'Total Align':<12} {'Prefix':<20}")
    click.echo("-" * 100)
    
    for idx, row in results.iterrows():
        click.echo(f"{idx:<6} {row['reference_org']:<40} {row['min_len']:<8} {row['mapped_contigs']:<8} {row['total_alignments']:<12} {row['prefix']:<20}")

@cli.command()
@click.argument('index', type=int)
@click.option('--browser', is_flag=True, help='Open in web browser')
@click.option('--output', type=click.Path(), help='Save HTML to file')
@click.option('--max-contigs', type=int, default=100, help='Maximum number of contigs to display')
@click.option('--min-length', type=int, default=0, help='Minimum alignment length')
def show(index, browser, output, max_contigs, min_length):
    """Show alignment details for specific index"""
    browser_obj = get_browser()
    if browser_obj.df is None:
        click.echo("❌ No data loaded. Use 'load' command first.")
        return
    
    record = browser_obj.get_record_by_index(index)
    if not record:
        click.echo(f"❌ Invalid index: {index}")
        return
    
    click.echo(f"\n📋 RECORD DETAILS - Index {index}")
    click.echo("=" * 50)
    click.echo(f"Taxon: {record['reference_org']}")
    click.echo(f"Min Length: {record['min_len']}")
    click.echo(f"Mapped Contigs: {record['mapped_contigs']}")
    click.echo(f"Total Alignments: {record['total_alignments']}")
    click.echo(f"Prefix: {record['prefix']}")
    click.echo(f"Log Path: {record['log_path']}")
    
    if record['mapped_contigs'] > 0:
        bam_path = Path(record['log_path']).parent / 'alignments.sorted.bam'
        if bam_path.exists():
            click.echo(f"\n🎨 Generating alignment plot...")
            fig = browser_obj.plot_alignments_plotly(
                bam_path, 
                title=record['reference_org'],
                max_contigs=max_contigs,
                min_length=min_length
            )
            
            if fig:
                if output:
                    fig.write_html(output)
                    click.echo(f"✅ Plot saved to: {output}")
                
                if browser or not output:
                    # Create temporary HTML file and open in browser
                    with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as f:
                        temp_file = f.name
                    fig.write_html(temp_file)
                    webbrowser.open(f'file://{temp_file}')
                    click.echo(f"✅ Plot opened in browser: {temp_file}")
            else:
                click.echo("❌ Failed to generate plot")
        else:
            click.echo("❌ BAM file not found")
    else:
        click.echo("⚠️  No mapped contigs to display")

@cli.command()
def interactive():
    """Start interactive browser session"""
    browser = get_browser()
    if browser.df is None:
        click.echo("❌ No data loaded. Use 'load' command first.")
        return
    
    click.echo("\n🎮 INTERACTIVE MINIMAP2 BROWSER")
    click.echo("=" * 50)
    click.echo("Commands:")
    click.echo("stats    - Show taxon statistics")
    click.echo("search   - Search for taxon")
    click.echo("show     - Show alignment details")
    click.echo("quit     - Exit browser")
    click.echo("=" * 50)
    
    while True:
        try:
            command = click.prompt("\nEnter command", type=str).strip().lower()
            
            if command == 'quit':
                save_browser()
                break
            elif command == 'stats':
                # ... (stats code remains the same)
                pass
            elif command == 'search':
                # ... (search code remains the same)
                pass
            elif command == 'show':
                index = click.prompt("Enter record index", type=int)
                record = browser.get_record_by_index(index)
                if record:
                    click.echo(f"\n📋 RECORD DETAILS - Index {index}")
                    click.echo("=" * 50)
                    click.echo(f"Taxon: {record['reference_org']}")
                    click.echo(f"Min Length: {record['min_len']}")
                    click.echo(f"Mapped Contigs: {record['mapped_contigs']}")
                    click.echo(f"Total Alignments: {record['total_alignments']}")
                    
                    if record['mapped_contigs'] > 0:
                        bam_path = Path(record['log_path']).parent / 'alignments.sorted.bam'
                        if bam_path.exists():
                            click.echo(f"\n🎨 Generating alignment plot...")
                            fig = browser.plot_alignments_plotly(
                                bam_path, 
                                title=record['reference_org'],
                                max_contigs=100
                            )
                            if fig:
                                # ✅ ИСПРАВЛЕНИЕ: Создаем осмысленное имя файла
                                safe_prefix = record['prefix'].replace('/', '_').replace('\\', '_').replace(' ', '_')
                                safe_taxon = record['reference_org'].replace('/', '_').replace('\\', '_').replace(' ', '_')
                                output_filename = f"{safe_prefix}_{safe_taxon}.html"
                                
                                # Сохраняем в текущую директорию
                                fig.write_html(output_filename)
                                webbrowser.open(f'file://{os.path.abspath(output_filename)}')
                                click.echo(f"✅ Plot saved and opened: {output_filename}")
                            else:
                                click.echo("❌ Failed to generate plot")
                        else:
                            click.echo("❌ BAM file not found")
                    else:
                        click.echo("⚠️  No mapped contigs to display")
                else:
                    click.echo(f"❌ Invalid index: {index}")
            else:
                click.echo("❌ Unknown command. Available: stats, search, show, quit")
                
        except KeyboardInterrupt:
            click.echo("\n👋 Goodbye!")
            save_browser()
            break
        except Exception as e:
            click.echo(f"❌ Error: {e}")

@cli.command()
def info():
    """Show current browser information"""
    browser = get_browser()
    if browser.df is None:
        click.echo("❌ No data loaded. Use 'load' command first.")
        return
    
    click.echo("\n📊 BROWSER INFORMATION")
    click.echo("=" * 30)
    click.echo(f"Total records: {len(browser.df)}")
    click.echo(f"Base pattern: {browser.base_pattern}")
    click.echo(f"Unique taxa: {browser.df['reference_org'].nunique()}")
    click.echo(f"Total mapped contigs: {browser.df['mapped_contigs'].sum()}")
    click.echo(f"Records with mappings: {(browser.df['mapped_contigs'] > 0).sum()}")

# Main execution
if __name__ == '__main__':
    cli()
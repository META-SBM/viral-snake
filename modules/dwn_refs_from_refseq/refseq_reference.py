import pandas as pd
from pathlib import Path
import subprocess
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

def find_and_combine_best_hits(prefix, min_len, preset, database):
    """
    Find all diamond hits files and combine best hits from multiple collections
    
    Parameters:
    prefix (str): Base path pattern
    min_len (int): Minimum contig length
    preset (str): Diamond preset (faster, sensitive, etc.)
    database (str): Database name (NR, SwissProt, etc.)
    """
    # Construct search pattern
    diamond_files  = list(Path(prefix).rglob(f"contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/hits_with_taxonomy.tsv"))
    
    print(f"Found {len(diamond_files)} collections with Diamond hits")
    
    all_best_hits = []
    
    for diamond_file in diamond_files:
        # Extract collection name from path (adjust index based on your directory structure)
        path_parts = diamond_file.parts
        collection = path_parts[-5]  # Adjust this index based on your directory structure
        
        print(f"Processing: {collection}")
        
        try:
            # Construct paths using the same prefix pattern
            stats_file = Path(prefix).rglob(f"contigs_formatted_minlen_{min_len}/contig_stats.tsv")
            
            # Read data
            stats = pd.read_csv(stats_file, sep='\t')
            diamond_hits = pd.read_csv(diamond_file, sep='\t')
            
            # Merge and get best hits
            merged_data = diamond_hits.merge(stats, left_on='qseqid', right_on='contig_id')
            best_hits = merged_data.loc[merged_data.groupby('qseqid')['bitscore'].idxmax()]
            
            # Add collection identifier
            best_hits['collection'] = collection
            
            all_best_hits.append(best_hits)
            
        except Exception as e:
            print(f"Error processing {collection}: {e}")
            continue
    
    if not all_best_hits:
        print("No data found!")
        return pd.DataFrame()
    
    # Combine all results
    combined_df = pd.concat(all_best_hits, ignore_index=True)
    
    print(f"\nCombined results:")
    print(f"Total rows: {len(combined_df)}")
    print(f"Collections: {combined_df['collection'].nunique()}")
    print(f"Collections found: {combined_df['collection'].unique()}")
    
    return combined_df

def filter_viral_hits(combined_df, min_identity=80, taxonomy_domain='Viruses'):
    """
    Filter results for viral hits with high identity
    
    Parameters:
    combined_df (pd.DataFrame): Combined Diamond hits dataframe
    min_identity (float): Minimum percent identity threshold
    taxonomy_domain (str): Taxonomy domain to filter (default: Viruses)
    """
    # Filter for specified domain
    domain_hits = combined_df[combined_df['Domain'] == taxonomy_domain]
    
    # Get best hit per contig and apply identity filter
    filtered_hits = domain_hits.loc[domain_hits.groupby('qseqid')['bitscore'].idxmax()]
    filtered_hits = filtered_hits[filtered_hits['pident'] > min_identity]
    
    print(f"Viral hits after filtering:")
    print(f"Total viral contigs: {len(filtered_hits)}")
    print(f"Unique species: {filtered_hits['Species'].nunique()}")
    
    return filtered_hits

def download_refseq_genomes(species_list, output_dir="refseq_genomes"):
    """
    Download reference genomes for a list of species from NCBI datasets
    
    Parameters:
    species_list (list): List of species names to download
    output_dir (str): Output directory for downloaded genomes
    """
    os.makedirs(output_dir, exist_ok=True)
    
    downloaded_species = []
    failed_species = []
    
    for species in species_list:
        print(f"Downloading genomes for: {species}")
        
        # Create safe filename
        safe_filename = species.replace(" ", "_").replace("/", "_").replace("\\", "_")
        output_zip = os.path.join(output_dir, f"{safe_filename}.zip")
        
        # Skip if already downloaded and file is not empty
        if os.path.exists(output_zip) and os.path.getsize(output_zip) > 1000:
            print(f"✓ Already downloaded, skipping: {species}")
            downloaded_species.append(species)
            continue
        
        try:
            # Download using datasets command
            cmd = [
                'datasets', 'download', 'virus', 'genome', 'taxon',
                species, 
                '--complete-only',
                '--filename', output_zip
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                # Check if archive is not empty
                if os.path.exists(output_zip) and os.path.getsize(output_zip) > 1000:
                    print(f"✓ Successfully downloaded: {species}")
                    downloaded_species.append(species)
                    
                    # Extract archive
                    extract_dir = os.path.join(output_dir, safe_filename)
                    os.makedirs(extract_dir, exist_ok=True)
                    
                    subprocess.run(['unzip', '-q', '-o', output_zip, '-d', extract_dir], check=True)
                    print(f"✓ Extracted to: {extract_dir}")
                else:
                    print(f"✗ Empty archive for: {species}")
                    failed_species.append(species)
                    if os.path.exists(output_zip):
                        os.remove(output_zip)
            else:
                print(f"✗ Download error for {species}: {result.stderr}")
                failed_species.append(species)
                
        except subprocess.TimeoutExpired:
            print(f"✗ Timeout for: {species}")
            failed_species.append(species)
        except Exception as e:
            print(f"✗ Error for {species}: {e}")
            failed_species.append(species)
    
    return downloaded_species, failed_species

def download_refseq_from_diamond_results(diamond_results_file, output_dir="refseq_genomes", 
                                       min_identity=80, taxonomy_domain='Viruses'):
    """
    Main function to download genomes based on Diamond results
    
    Parameters:
    diamond_results_file (str): Path to Diamond results file
    output_dir (str): Output directory for reference genomes
    min_identity (float): Minimum percent identity for filtering
    taxonomy_domain (str): Taxonomy domain to filter
    """
    # Read Diamond results
    diamond_df = pd.read_csv(diamond_results_file, sep='\t')
    
    # Filter for viral hits
    viral_hits = filter_viral_hits(diamond_df, min_identity, taxonomy_domain)
    
    # Extract unique species
    species_column = None
    for col in ['Species', 'species', 'staxids', 'scientific_name']:
        if col in viral_hits.columns:
            species_column = col
            break
    
    if species_column is None:
        print("No species column found in the data!")
        return
    
    unique_species = viral_hits[species_column].dropna().unique()
    
    print(f"Found {len(unique_species)} unique species for download")
    
    # Download genomes
    downloaded, failed = download_refseq_genomes(unique_species, output_dir)
    
    # Print report
    print(f"\n=== DOWNLOAD REPORT ===")
    print(f"Successfully downloaded: {len(downloaded)}")
    print(f"Failed to download: {len(failed)}")
    
    if failed:
        print("\nFailed species:")
        for species in failed:
            print(f" - {species}")

def extract_protein_sequences_from_annotation(ref_fasta, protein_type_patterns):
    """
    Extract protein sequences from FASTA file based on name patterns
    and save in the same directory with genomic_ prefix
    
    Args:
        ref_fasta (str): path to FASTA file with reference sequences
        protein_type_patterns (list): list of tuples (pattern, protein_type)
    
    Returns:
        dict: dictionary with paths to created files {protein_type: file_path}
    """
    
    # Output directory - same as ref_fasta location
    output_dir = Path(ref_fasta).parent
    
    # Dictionary to store sequences by protein types
    protein_sequences = {protein_type: [] for _, protein_type in protein_type_patterns}
    protein_sequences['other'] = []  # for sequences not matching any pattern
    
    # Read FASTA file
    try:
        sequences = list(SeqIO.parse(ref_fasta, "fasta"))
        print(f"  Read {len(sequences)} sequences")
    except Exception as e:
        print(f"  Error reading file {ref_fasta}: {e}")
        return {}
    
    # Compile regex patterns for performance
    compiled_patterns = [(re.compile(pattern, re.IGNORECASE), protein_type) 
                        for pattern, protein_type in protein_type_patterns]
    
    # Process each sequence
    for record in sequences:
        sequence_found = False
        description = record.description
        
        # Check all patterns
        for pattern, protein_type in compiled_patterns:
            if pattern.search(description):
                # Create new record
                new_record = SeqRecord(
                    seq=record.seq,
                    id=record.id,
                    description=record.description  # keep original description
                )
                protein_sequences[protein_type].append(new_record)
                sequence_found = True
                break
        
        # If sequence doesn't match any pattern
        if not sequence_found:
            protein_sequences['other'].append(record)
    
    # Save sequences to separate files
    output_files = {}
    
    for protein_type, sequences in protein_sequences.items():
        if sequences:  # save only if there are sequences
            # Replace spaces with underscores in protein type name
            protein_type_safe = protein_type.replace(' ', '_')
            output_filename = f"genomic_{protein_type_safe}.fasta"
            output_path = output_dir / output_filename
            
            try:
                SeqIO.write(sequences, output_path, "fasta")
                output_files[protein_type] = output_path
                print(f"  Created file {output_filename} with {len(sequences)} sequences")
            except Exception as e:
                print(f"  Error writing file {output_path}: {e}")
    
    return output_files

def find_all_refseq_taxa(refseq_base_dir):
    """
    Find all taxa in refseq_reference directory
    
    Args:
        refseq_base_dir (str): base directory of refseq_reference
    
    Returns:
        list: list of paths to genomic.fna files for each taxon
    """
    refseq_path = Path(refseq_base_dir)
    genomic_files = []
    
    if not refseq_path.exists():
        print(f"Directory {refseq_base_dir} does not exist!")
        return genomic_files
    
    # Find all genomic.fna files in subdirectories
    for genomic_file in refseq_path.rglob("*/ncbi_dataset/data/genomic.fna"):
        genomic_files.append(genomic_file)
    
    return genomic_files

def process_single_taxon(genomic_file, protein_patterns):
    """
    Process single taxon
    
    Returns:
        tuple: (taxon_name, success_bool, result_files)
    """
    taxon_name = genomic_file.parent.parent.parent.name
    print(f"Processing taxon: {taxon_name}")
    
    try:
        # Extract proteins - files are saved in the same directory as genomic.fna
        result_files = extract_protein_sequences_from_annotation(
            ref_fasta=str(genomic_file),
            protein_type_patterns=protein_patterns
        )
        
        return (taxon_name, True, result_files)
        
    except Exception as e:
        print(f"Error processing taxon {taxon_name}: {e}")
        return (taxon_name, False, {})

def process_all_refseq_taxa(refseq_base_dir, protein_patterns, max_workers=4):
    """
    Process all taxa in refseq_reference
    
    Args:
        refseq_base_dir (str): path to refseq_reference directory
        protein_patterns (list): list of patterns for protein classification
        max_workers (int): number of parallel processes
    
    Returns:
        dict: processing statistics
    """
    
    # Find all genomic.fna files
    print("Searching for taxa in refseq_reference...")
    genomic_files = find_all_refseq_taxa(refseq_base_dir)
    
    if not genomic_files:
        print("No genomic.fna files found!")
        return {}
    
    print(f"Found taxa: {len(genomic_files)}")
    
    # Statistics
    stats = {
        'total': len(genomic_files),
        'success': 0,
        'failed': 0,
        'results': {}
    }
    
    # Process taxa
    print("Starting taxa processing...")
    
    # Sequential processing
    for genomic_file in tqdm(genomic_files, desc="Processing taxa"):
        taxon_name, success, result_files = process_single_taxon(
            genomic_file, protein_patterns
        )
        
        stats['results'][taxon_name] = {
            'success': success,
            'result_files': result_files
        }
        
        if success:
            stats['success'] += 1
        else:
            stats['failed'] += 1
    
    return stats

def print_protein_extraction_statistics(stats):
    """Print protein extraction statistics"""
    print("\n" + "="*50)
    print("PROTEIN EXTRACTION STATISTICS")
    print("="*50)
    print(f"Total taxa: {stats['total']}")
    print(f"Successfully processed: {stats['success']}")
    print(f"Failed to process: {stats['failed']}")
    
    # Count total extracted files by protein types
    protein_counts = {}
    for taxon, result in stats['results'].items():
        if result['success']:
            for protein_type, file_path in result['result_files'].items():
                if protein_type not in protein_counts:
                    protein_counts[protein_type] = 0
                # Count sequences in file
                try:
                    count = len(list(SeqIO.parse(file_path, "fasta")))
                    protein_counts[protein_type] += count
                except:
                    pass
    
    print("\nExtracted sequences by protein types:")
    for protein_type, count in sorted(protein_counts.items()):
        print(f"  {protein_type}: {count} sequences")

def check_results_for_taxon(taxon_path):
    """Check results for specific taxon"""
    taxon_dir = Path(taxon_path)
    data_dir = taxon_dir / "ncbi_dataset" / "data"
    
    if not data_dir.exists():
        print(f"Directory not found: {data_dir}")
        return
    
    print(f"\nFiles in {data_dir}:")
    for file in sorted(data_dir.iterdir()):
        if file.is_file():
            size = file.stat().st_size
            print(f"  {file.name} ({size} bytes)")

# Complete workflow function
def run_complete_workflow(config):
    """
    Run complete workflow: Diamond analysis → Genome download → Protein extraction
    
    Parameters:
    config (dict): Configuration dictionary with all parameters
    """
    
    print("="*60)
    print("STARTING COMPplete VIRAL ANALYSIS WORKFLOW")
    print("="*60)
    
    # Step 1: Diamond analysis
    print("\n1. DIAMOND ANALYSIS")
    print("-" * 30)
    
    combined_results = find_and_combine_best_hits(
        prefix=config['diamond_prefix'],
        min_len=config['min_len'],
        preset=config['diamond_preset'],
        database=config['database']
    )
    
    if combined_results.empty:
        print("No Diamond results found! Exiting workflow.")
        return
    
    # Save combined results
    combined_results.to_csv(config['diamond_output'], sep='\t', index=False)
    print(f"Combined results saved to: {config['diamond_output']}")
    
    # Step 2: Download reference genomes
    print("\n2. REFERENCE GENOME DOWNLOAD")
    print("-" * 30)
    
    download_refseq_from_diamond_results(
        diamond_results_file=config['diamond_output'],
        output_dir=config['refseq_dir'],
        min_identity=config['min_identity'],
        taxonomy_domain=config['taxonomy_domain']
    )
    
    # Step 3: Protein sequence extraction
    print("\n3. PROTEIN SEQUENCE EXTRACTION")
    print("-" * 30)
    
    # Protein patterns for classification
    protein_type_patterns = [
        (r'\b(RNA-dependent RNA-polymerase|RdRp protein|polymerase|RNA polymerase)\b', 'polymerase'),
        (r'\b(nucleocapsid protein|nucleocapsid|capsid)\b', 'nucleocapsid_protein'),
        (r'\b(tail protein|tail)\b', 'tail_protein'),
        (r'\b(L protein|large protein)\b', 'L_protein'),
        (r'\b(NS3-like|NS3|nonstructural protein 3)\b', 'NS3-like'),
        (r'\b(Toprim-like|Toprim|primase)\b', 'Toprim-like'),
        (r'\b(glycoprotein|glycoprotein precursor)\b', 'glycoprotein'),
        (r'\b(NS5|nonstructural protein 5)\b', 'NS5-like'),
        (r'\b(helicase)\b', 'helicase'),
        (r'\b(methyltransferase)\b', 'methyltransferase'),
        (r'\b(protease)\b', 'protease'),
        (r'\b(NS2|nonstructural protein 2)\b', 'NS2-like'),
        (r'\b(NS1|nonstructural protein 1)\b', 'NS1-like'),
        (r'\b(VP1|VP2|VP3|VP4|viral protein)\b', 'viral_protein'),
        (r'\b(envelope protein|envelope)\b', 'envelope_protein'),
        (r'\b(matrix protein|matrix)\b', 'matrix_protein'),
        (r'\b(segment S)\b', 'segment_S'),
        (r'\b(segment L)\b', 'segment_L'),
        (r'\b(segment M)\b', 'segment_M'),
        (r'\b(segment 4)\b', 'segment_4'),
        (r'\b(segment 3)\b', 'segment_3'),
        (r'\b(segment 2)\b', 'segment_2'),
        (r'\b(segment 1)\b', 'segment_1'),
        (r'.+', 'other')
    ]
    
    stats = process_all_refseq_taxa(
        refseq_base_dir=config['refseq_dir'],
        protein_patterns=protein_type_patterns,
        max_workers=config.get('max_workers', 4)
    )
    
    # Print protein extraction statistics
    print_protein_extraction_statistics(stats)
    
    print("\n" + "="*60)
    print("WORKFLOW COMPLETED SUCCESSFULLY")
    print("="*60)

# Example configuration
if __name__ == "__main__":
    config = {
        'diamond_prefix': '/mnt/mgx/DATASETS/INTERNAL/VIROME/RUN3/co_assembly/megahit',
        'min_len': 800,
        'diamond_preset': 'faster',
        'database': 'NR',
        'min_identity': 80,
        'taxonomy_domain': 'Viruses',
        'diamond_output': '/mnt/mgx/RUNS/kuzmichenko_pa/third_run/diamond_combined_results.tsv',
        'refseq_dir': '/mnt/mgx/DATASETS/INTERNAL/VIROME/RUN3/refseq_reference',
        'max_workers': 4
    }
    
    run_complete_workflow(config)
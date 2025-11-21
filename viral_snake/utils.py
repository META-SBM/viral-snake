"""Utility functions for viral-snake pipeline"""

import os
import glob
from pathlib import Path
import pandas as pd


# ============================================================================
# Sample Discovery
# ============================================================================

def get_samples(reads_dir, r1_pattern="_R1.fastq.gz", r2_pattern="_R2.fastq.gz"):
    """
    Extract sample names from paired-end read files.
    
    Args:
        reads_dir: Directory containing read files
        r1_pattern: Pattern to identify R1 files
        r2_pattern: Pattern to identify R2 files (for validation)
    
    Returns:
        list: Sorted sample names
        
    Raises:
        FileNotFoundError: If R1/R2 pairs don't match
    """
    reads_path = Path(reads_dir)
    
    # Get R1 files
    r1_files = sorted(reads_path.glob(f"*{r1_pattern}"))
    samples = [f.name.replace(r1_pattern, "") for f in r1_files]
    
    # Validate R2 files exist
    missing_r2 = []
    for sample in samples:
        r2_file = reads_path / f"{sample}{r2_pattern}"
        if not r2_file.exists():
            missing_r2.append(str(r2_file))
    
    if missing_r2:
        raise FileNotFoundError(
            f"Missing R2 files for paired samples:\n" + "\n".join(missing_r2)
        )
    
    return samples


def get_samples_from_manifest(manifest_file, sample_col="sample_id"):
    """
    Load sample names from a TSV/CSV manifest file.
    
    Args:
        manifest_file: Path to manifest file
        sample_col: Column name containing sample IDs
    
    Returns:
        list: Sample names
    """
    df = pd.read_csv(manifest_file, sep="\t" if manifest_file.endswith('.tsv') else ",")
    return df[sample_col].tolist()


# ============================================================================
# Reference Discovery
# ============================================================================

def discover_references(refseq_dir, pattern='*/ncbi_dataset/data/genomic.fna'):
    """
    Discover reference organisms from refseq directory structure.
    
    Args:
        refseq_dir: Base directory containing reference genomes
        pattern: Glob pattern relative to refseq_dir
    
    Returns:
        list: Sorted reference organism identifiers
    """
    refs = glob.glob(os.path.join(refseq_dir, pattern))
    # Extract organism name from path (4th from end: /org_name/ncbi_dataset/data/genomic.fna)
    ref_names = [Path(r).parts[-4] for r in refs]
    return sorted(ref_names)


# ============================================================================
# File Validation
# ============================================================================

def validate_paired_reads(reads_dir, samples, r1_suffix="_R1.fastq.gz", r2_suffix="_R2.fastq.gz"):
    """
    Validate that both R1 and R2 files exist for all samples.
    
    Args:
        reads_dir: Directory containing reads
        samples: List of sample names
        r1_suffix: R1 file suffix
        r2_suffix: R2 file suffix
    
    Returns:
        bool: True if all files exist
        
    Raises:
        FileNotFoundError: If any files are missing
    """
    reads_path = Path(reads_dir)
    missing = []
    
    for sample in samples:
        r1 = reads_path / f"{sample}{r1_suffix}"
        r2 = reads_path / f"{sample}{r2_suffix}"
        
        if not r1.exists():
            missing.append(str(r1))
        if not r2.exists():
            missing.append(str(r2))
    
    if missing:
        raise FileNotFoundError(
            f"Missing {len(missing)} read file(s):\n" + "\n".join(missing)
        )
    
    return True

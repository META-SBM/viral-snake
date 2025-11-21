"""Collection management for co-assembly"""

import yaml
from pathlib import Path


def load_collections_from_dir(collection_dir):
    """Load all collection YAML files from a directory
    
    Args:
        collection_dir: Path to directory containing collection YAML files
        
    Returns:
        dict: Collections indexed by name
    """
    collections = {}
    collection_path = Path(collection_dir)
    
    if not collection_path.exists():
        return collections
    
    for yaml_file in collection_path.glob("*.yaml"):
        with open(yaml_file) as f:
            coll = yaml.safe_load(f)
            if 'name' in coll:
                collections[coll['name']] = coll
            else:
                print(f"Warning: Collection in {yaml_file} has no 'name' field")
    
    return collections


def load_collections(config_file):
    """Load collection definitions from YAML config
    
    Returns:
        dict: Collections with format:
            {
                'collection_name': {
                    'samples': ['sample1', 'sample2'],
                    'qc_filter': 'raw__cutadapt_no_mgi_min_len_90'
                }
            }
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    return config.get('collections', {})


def validate_collections(collections, available_samples):
    """Validate that all samples in collections exist"""
    errors = []
    for coll_name, coll_data in collections.items():
        samples = coll_data.get('samples', [])
        for sample in samples:
            if sample not in available_samples:
                errors.append(f"Collection '{coll_name}': sample '{sample}' not found")
    
    if errors:
        raise ValueError("\n".join(errors))
    return True
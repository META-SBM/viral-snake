if __name__ == "__main__":
    
    import argparse
    import pandas as pd
    import os
    from subprocess import run, Popen, PIPE
    from time import time

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-input', '--input_path', 
                        help='input directory path', 
                        required=True)
    parser.add_argument('-output', '--output_path', 
                        help='output directory path', 
                        required=True)
    parser.add_argument('-coef', 
                        help='trashold for target is bitscore * coef', 
                        type=float, 
                        required=True)
    parser.add_argument('-taxdump', '--taxdump_path',
                        help='path to taxdump directory',
                        required=True)

    args = vars(parser.parse_args())

    i_path = args['input_path']
    o_path = args['output_path']
    coef = args['coef']
    taxdump = args['taxdump_path']

    script_start = time()
    print(f"[DEBUG] Script started")
    
    read_start = time()
    diam_res = pd.read_csv(i_path, sep='\t').dropna(subset=['staxids'])
    print(f"[DEBUG] Read input file: {time() - read_start:.2f}s ({len(diam_res)} rows, {diam_res['qseqid'].nunique()} unique queries)")
    
    # Prepare all data first
    prep_start = time()
    query_ids = []
    lca_inputs = []
    tophit_inputs = []
    
    # Prepare all data first (VECTORIZED - FASTEST)
    prep_start = time()

    # Calculate max bitscore per query using transform (broadcast back to original df)
    diam_res['max_bitscore'] = diam_res.groupby('qseqid')['bitscore'].transform('max')

    # Filter all rows at once
    filtered_df = diam_res[diam_res['bitscore'] >= diam_res['max_bitscore'] * coef].copy()

    query_ids = []
    lca_inputs = []
    tophit_inputs = []

    # Single groupby iteration - ensures order consistency
    for qseqid, group in filtered_df.groupby('qseqid'):
        query_ids.append(qseqid)
        
        # Get row with maximum bitscore
        top_hit_row = group.loc[group['bitscore'].idxmax()]
        tophit_inputs.append(str(top_hit_row['staxids']))
        
        # Collect all taxonomy IDs from all hits for LCA calculation
        # Example: hits have staxids ["562;561", "562", "1613;561"]
        #          Result: ["562", "561", "562", "1613", "561"]
        taxids = []
        for staxids in group['staxids'].dropna().values:
            # Split "562;561" → ["562", "561"] and add to list
            taxids.extend(str(staxids).split(';'))
        
        # Format for taxonkit: ["562", "561", "543"] → "562;561;543"
        lca_inputs.append(';'.join(taxids))

    print(f"[DEBUG] Data preparation: {time() - prep_start:.2f}s")
    
    # BATCH PROCESS ALL LCA QUERIES AT ONCE
    lca_start = time()
    lca_batch_input = '\n'.join(lca_inputs)
    result = run(f"taxonkit lca --data-dir {taxdump} -s ';' | cut -f2 | taxonkit reformat2 --data-dir {taxdump} -I 1 -r 'unknown' -f '{{domain|superkingdom}};{{phylum}};{{class}};{{order}};{{family}};{{genus}};{{species}};{{subspecies|strain|no rank}}'", 
                shell=True, input=lca_batch_input, capture_output=True, text=True)
    lca_results = result.stdout.strip().split('\n')
    print(f"[DEBUG] Batch LCA processing: {time() - lca_start:.2f}s ({len(lca_inputs)} queries)")
    
    # BATCH PROCESS ALL TOP HITS AT ONCE
    tophit_start = time()
    tophit_batch_input = '\n'.join(tophit_inputs)
    result = run(f"taxonkit lca --data-dir {taxdump} -s ';' | cut -f2 | taxonkit reformat2 --data-dir {taxdump} -I 1 -r 'unknown' -f '{{domain|superkingdom}};{{phylum}};{{class}};{{order}};{{family}};{{genus}};{{species}};{{subspecies|strain|no rank}}'", 
                shell=True, input=tophit_batch_input, capture_output=True, text=True)
    tophit_results = result.stdout.strip().split('\n')
    print(f"[DEBUG] Batch top hit processing: {time() - tophit_start:.2f}s ({len(tophit_inputs)} queries)")
    
    # Assemble results: each query becomes one row
    assembly_start = time()
    rows = []  # List to hold all rows

    for idx, qseqid in enumerate(query_ids):
        # Parse LCA result: "taxid\ttaxonomy"
        lca_parts = lca_results[idx].split('\t')
        lca_taxid = lca_parts[0]
        lca_taxonomy = lca_parts[1] if len(lca_parts) > 1 else 'unknown'
        taxonomy_parts = lca_taxonomy.split(';')
        
        # Parse top hit result: "taxid\tfull_lineage"
        tophit_parts = tophit_results[idx].split('\t')
        tophit_taxid = tophit_parts[0]
        tophit_taxonomy = tophit_parts[1] if len(tophit_parts) > 1 else 'unknown'
        
        row = {
            'qseqid': qseqid,
            'LCA_taxid': lca_taxid,
            'Domain': taxonomy_parts[0] if len(taxonomy_parts) > 0 else 'unknown',
            'Phylum': taxonomy_parts[1] if len(taxonomy_parts) > 1 else 'unknown',
            'Class': taxonomy_parts[2] if len(taxonomy_parts) > 2 else 'unknown',
            'Order': taxonomy_parts[3] if len(taxonomy_parts) > 3 else 'unknown',
            'Family': taxonomy_parts[4] if len(taxonomy_parts) > 4 else 'unknown',
            'Genus': taxonomy_parts[5] if len(taxonomy_parts) > 5 else 'unknown',
            'Species': taxonomy_parts[6] if len(taxonomy_parts) > 6 else 'unknown',
            'Subsp': taxonomy_parts[7] if len(taxonomy_parts) > 7 else 'unknown',
            'TopHit_taxid': tophit_taxid,
            'TopHit': tophit_taxonomy 
        }
        
        rows.append(row)

    print(f"[DEBUG] Result assembly: {time() - assembly_start:.2f}s")

    # Convert list of dicts to DataFrame
    output_start = time()
    LCA_res = pd.DataFrame(rows)
    LCA_res.set_index('qseqid', inplace=True)
    LCA_res.to_csv(o_path, sep='\t')
    print(f"[DEBUG] Output written: {time() - output_start:.2f}s")
    
    total_time = time() - script_start
    print(f"[DEBUG] Total runtime: {total_time:.2f}s")
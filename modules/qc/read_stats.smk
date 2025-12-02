rule count_reads:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        stats = "qc/read_stats/{qc_filter}/{sample}_read_counts.tsv"
    threads: 4
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        # Count reads in both files
        seqkit stats -j {threads} -T {input.r1} {input.r2} > {output.stats}
        """

def get_collection_read_stats(wildcards):
    """Get read stats files for all samples in a collection"""
    import yaml
    
    # Load collection metadata
    collection_file = f"config/sample_collections/{wildcards.collection}.yaml"
    with open(collection_file) as f:
        collection = yaml.safe_load(f)
    
    qc_filter = collection['qc_filter']
    samples = collection['samples']
    
    # Return list of stats files for these samples
    return expand(
        "qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
        qc_filter=qc_filter,
        sample=samples
    )


rule collection_read_stats:
    """Aggregate read statistics for a sample collection"""
    input:
        collection_meta = "config/sample_collections/{collection}.yaml",
        stats_files = get_collection_read_stats
    output:
        combined_stats = "qc/read_stats/collections/{collection}_read_stats.tsv"
    run:
        import yaml
        import pandas as pd
        from pathlib import Path
        
        # Load collection metadata
        with open(input.collection_meta) as f:
            collection = yaml.safe_load(f)
        
        qc_filter = collection['qc_filter']
        samples = collection['samples']
        
        # Read all stats files and combine
        dfs = []
        for stats_file in input.stats_files:
            # Read seqkit stats output
            df = pd.read_csv(stats_file, sep='\t')
            
            # Filter to R1 only
            df_r1 = df[df['file'].str.contains('_R1.fastq.gz')].copy()
            
            # Extract sample name from file path
            # e.g., qc/read_stats/raw/sample123_read_counts.tsv -> sample123
            sample_name = Path(stats_file).stem.replace('_read_counts', '')
            
            # Add columns at the beginning
            df_r1.insert(0, 'qc_filter', qc_filter)
            df_r1.insert(1, 'sample', sample_name)
            
            # Clean up file column to just show filename
            df_r1['file'] = df_r1['file'].apply(lambda x: Path(x).name)
            
            dfs.append(df_r1)
        
        # Combine all dataframes
        combined = pd.concat(dfs, ignore_index=True)
        
        # Sort by sample name
        combined = combined.sort_values('sample')
        
        # Write output
        combined.to_csv(output.combined_stats, sep='\t', index=False)
        
        # Print summary
        total_reads = combined['num_seqs'].sum()
        total_bases = combined['sum_len'].sum()
        
        print(f"\n{'='*50}")
        print(f"Collection Read Statistics")
        print(f"{'='*50}")
        print(f"Collection:   {wildcards.collection}")
        print(f"QC Filter:    {qc_filter}")
        print(f"Samples:      {len(samples)}")
        print(f"Total reads:  {total_reads:,}")
        print(f"Total bases:  {total_bases:,}")
        print(f"Avg reads:    {total_reads/len(samples):,.0f} per sample")
        print(f"Output:       {output.combined_stats}")
        print(f"{'='*50}\n")
rule kraken2_classify:
    """
    Classify reads taxonomically using Kraken2.
    Generates both detailed output and summary report.
    """
    input:
        r1 = "{fs_prefix}/{dataset}/reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "{fs_prefix}/{dataset}/reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        report = "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.kraken2.report",
        output = "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.kraken2.output"
    params:
        db = DATABASES['kraken2']
    benchmark:
        "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.kraken2.benchmark.txt"
    log:
        "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.kraken2.log"
    threads: THREADS['kraken2']
    wildcard_constraints:
        dataset = "[^/]+",
        sample = "[^/]+",
        confidence = "[0-9.]+"
    conda:
        "../../envs/kraken2.yaml"
    shell:
        """
        echo "=== Running Kraken2 classification ===" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        echo "Sample: {wildcards.sample}" >> {log}
        echo "Confidence: {wildcards.confidence}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "" >> {log}
        
        k2 classify --db {params.db} \
            --threads {threads} \
            --paired \
            --confidence {wildcards.confidence} \
            --memory-mapping \
            --report {output.report} \
            --output {output.output} \
            --log {log} \
            {input.r1} {input.r2}
        """


rule bracken_abundance:
    """
    Estimate species-level abundance using Bracken.
    Re-estimates abundances from Kraken2 results using Bayesian methods.
    """
    input:
        report = "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.kraken2.report"
    output:
        bracken = "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.bracken",
        report = "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.bracken.report"
    params:
        db = DATABASES['kraken2'],
        read_len = 150,  # CHANGE THIS to match your actual read length
        level = "S",     # S=Species, G=Genus, F=Family, etc.
        threshold = 1    # Minimum number of reads required
    threads: 1
    wildcard_constraints:
        dataset = "[^/]+",
        sample = "[^/]+",
        confidence = "[0-9.]+"
    log:
        "{fs_prefix}/{dataset}/kraken2/{confidence}/{qc_filter}/{sample}.bracken.log"
    conda:
        "../../envs/kraken2.yaml"
    shell:
        """
        echo "Running Bracken abundance estimation" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        echo "Sample: {wildcards.sample}" >> {log}
        
        bracken -d {params.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold} \
            2>> {log}
        """


def get_bracken_table_inputs(wildcards):
    """
    Read meta.yaml and reconstruct bracken file paths for a feature table.
    
    Args:
        wildcards: Must contain fs_prefix, dataset, and feature_table_id
        
    Returns:
        list: Paths to bracken files for all samples in the feature table
    """
    import yaml
    
    meta_file = f"{wildcards.fs_prefix}/{wildcards.dataset}/feature_tables/bracken-{wildcards.feature_table_id}/meta.yaml"
    
    with open(meta_file, 'r') as f:
        meta = yaml.safe_load(f)
    
    # Reconstruct paths with fs_prefix/dataset
    confidence = meta.get('confidence', '0.5')
    qc_filter = meta['qc_filter']
    samples = meta['samples']
    
    return [
        f"{wildcards.fs_prefix}/{wildcards.dataset}/kraken2/{confidence}/{qc_filter}/{sample}.bracken"
        for sample in samples
    ]


def get_abundance_metric(wildcards):
    """
    Read meta.yaml to get abundance metric for feature table.
    
    Args:
        wildcards: Must contain fs_prefix, dataset, and feature_table_id
        
    Returns:
        str: Abundance metric name (e.g., 'new_est_reads', 'fraction_total_reads')
    """
    import yaml
    
    meta_file = f"{wildcards.fs_prefix}/{wildcards.dataset}/feature_tables/bracken-{wildcards.feature_table_id}/meta.yaml"
    
    with open(meta_file, 'r') as f:
        meta = yaml.safe_load(f)
    
    return meta.get('abundance_metric', 'new_est_reads')


rule build_bracken_abundance_table:
    """
    Build abundance matrix from multiple Bracken outputs.
    Creates a samples × taxa table for downstream analysis.
    """
    input:
        files = get_bracken_table_inputs
    output:
        table = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/abundance_table.tsv"
    params:
        metric = get_abundance_metric,
        input_args = lambda wildcards, input: ' '.join([f'-i {f}' for f in input.files]),
        script = "~/MGX/viral-snake/scripts/build_abundance_table.py"
    wildcard_constraints:
        dataset = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/build_table.log"
    shell:
        """
        echo "Building abundance table" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        echo "Feature table ID: {wildcards.feature_table_id}" >> {log}
        echo "Abundance metric: {params.metric}" >> {log}
        
        python {params.script} \
            {params.input_args} \
            --abundance-metric {params.metric} \
            --output {output.table} \
            2>> {log}
        """


rule create_taxonomy_table:
    """
    Create taxonomy reference table from abundance table tax IDs.
    Uses taxonkit to resolve tax IDs to full taxonomic lineages.
    """
    input:
        abundance = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/abundance_table.tsv"
    output:
        taxonomy = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/taxonomy_table.tsv"
    params:
        taxdump = DATABASES['taxdump']
    wildcard_constraints:
        dataset = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/taxonomy.log"
    conda:
        "../../envs/taxonkit.yaml"
    shell:
        """
        echo "Creating taxonomy table" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        
        # Extract tax IDs from first column (skip header, remove tax_id__ prefix)
        tail -n +2 {input.abundance} | cut -f1 | sed 's/tax_id__//' > {output.taxonomy}.tmp.ids
        
        echo "Extracted $(wc -l < {output.taxonomy}.tmp.ids) unique tax IDs" >> {log}
        
        # Get lineage and reformat with standard ranks
        cat {output.taxonomy}.tmp.ids | \
            taxonkit lineage --data-dir {params.taxdump} 2>> {log} | \
            taxonkit reformat2 \
                --data-dir {params.taxdump} \
                -I 1 \
                -f "{{domain|acellular root|superkingdom}};{{phylum}};{{class}};{{order}};{{family}};{{genus}};{{species}}" \
                -r "NA" \
                2>> {log} > {output.taxonomy}.tmp.reformatted
        
        # Extract taxid (col 1) and reformatted lineage (col 3)
        # Split semicolon-separated lineage into separate columns
        cut -f1,3 {output.taxonomy}.tmp.reformatted | \
            awk 'BEGIN {{FS="\t"; OFS="\t"}} 
                 {{
                   split($2, ranks, ";");
                   print "tax_id__"$1, ranks[1], ranks[2], ranks[3], ranks[4], ranks[5], ranks[6], ranks[7]
                 }}' > {output.taxonomy}.tmp.prefixed
        
        # Add header
        echo -e "tax_id\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" | \
            cat - {output.taxonomy}.tmp.prefixed > {output.taxonomy}
        
        # Clean up
        rm {output.taxonomy}.tmp.*
        
        echo "Taxonomy table complete" >> {log}
        """
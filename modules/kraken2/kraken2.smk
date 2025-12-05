from pathlib import Path

# Define script directory
KRAKEN2_MODULE_DIR = Path(workflow.basedir) / "modules/kraken2"
# Load heatmap presets
with open(KRAKEN2_MODULE_DIR / "heatmap_presets.yaml") as f:
    HEATMAP_PRESETS = yaml.safe_load(f)['heatmaps']

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
        script = KRAKEN2_MODULE_DIR / "build_abundance_table.py" 
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

def build_heatmap_command(wildcards):
    """
    Build heatmap command from preset configuration.
    
    Args:
        wildcards: Must contain heatmap_preset
        
    Returns:
        str: Command-line arguments for create_heatmap.R
    """
    preset = HEATMAP_PRESETS[wildcards.heatmap_preset]
    
    cmd_parts = [
        f"--prevalence {preset.get('prevalence', 0.1)}",
        f"--detection {preset.get('detection', 0)}",
    ]
    
    if 'domain' in preset:
        cmd_parts.append(f"--domain '{preset['domain']}'")
    
    return " ".join(cmd_parts)

def get_phyloseq_metadata_arg(wildcards):
    """Check if metadata file exists and return argument."""
    metadata_file = f"{wildcards.fs_prefix}/{wildcards.dataset}/feature_tables/bracken-{wildcards.feature_table_id}/metadata.tsv"
    if Path(metadata_file).exists():
        return f"-m {metadata_file}"
    return ""


rule create_phyloseq:
    """
    Create phyloseq object from abundance and taxonomy tables.
    Optionally includes sample metadata if available.
    """
    input:
        abundance = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/abundance_table.tsv",
        taxonomy = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/taxonomy_table.tsv"
    output:
        phyloseq = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/phyloseq.rds"
    params:
        script = KRAKEN2_MODULE_DIR / "create_phyloseq.R",
        metadata_arg = get_phyloseq_metadata_arg
    wildcard_constraints:
        dataset = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/phyloseq.log"
    conda:
        "../../envs/phyloseq.yaml"
    shell:
        """
        Rscript {params.script} \
            -a {input.abundance} \
            -t {input.taxonomy} \
            {params.metadata_arg} \
            -o {output.phyloseq} \
            2>&1 | tee {log}
        """

rule create_abundance_heatmap:
    """
    Create filtered abundance heatmap from phyloseq object.
    Filters taxa by prevalence and creates log-transformed heatmap using presets.
    """
    input:
        phyloseq = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/phyloseq.rds"
    output:
        heatmap = "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/heatmap_{heatmap_preset}.pdf"
    params:
        script = KRAKEN2_MODULE_DIR / "plot_heatmap.R",
        heatmap_cmd = build_heatmap_command
    wildcard_constraints:
        dataset = "[^/]+",
        heatmap_preset = "|".join(HEATMAP_PRESETS.keys())
    log:
        "{fs_prefix}/{dataset}/feature_tables/bracken-{feature_table_id}/heatmap_{heatmap_preset}.log"
    conda:
        "../../envs/phyloseq.yaml"
    shell:
        """
        Rscript {params.script} \
            -p {input.phyloseq} \
            -o {output.heatmap} \
            {params.heatmap_cmd} \
            2>&1 | tee {log}
        """

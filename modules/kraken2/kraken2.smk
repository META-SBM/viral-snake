rule kraken2_classify:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        report = "kraken2/{qc_filter}/{sample}.kraken2.report",
        output = "kraken2/{qc_filter}/{sample}.kraken2.output"
    params:
        db = DATABASES['kraken2'],
        confidence = 0.4
    log:
        "kraken2/{qc_filter}/{sample}.kraken2.log"
    threads: 16
    conda:
        "../../envs/kraken2.yaml"
    shell:
        """
        k2 classify --db {params.db} \
            --threads {threads} \
            --paired \
            --confidence {params.confidence} \
            --memory-mapping \
            --report {output.report} \
            --output {output.output} \
            --log {log} \
            {input.r1} {input.r2}
        """


rule bracken_abundance:
    input:
        report = "kraken2/{qc_filter}/{sample}.kraken2.report"
    output:
        bracken = "kraken2/{qc_filter}/{sample}.bracken",
        report = "kraken2/{qc_filter}/{sample}.bracken.report"
    params:
        db = DATABASES['kraken2'],
        read_len = 150,  # CHANGE THIS to match your actual read length
        level = "S",     # S=Species, G=Genus, F=Family, etc.
        threshold = 10   # Minimum number of reads required
    threads: 1
    conda:
        "../../envs/kraken2.yaml"
    shell:
        """
        bracken -d {params.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold}
        """


def get_bracken_table_inputs(wildcards):
    """Read meta.yaml and reconstruct bracken file paths."""
    meta_file = f"feature_tables/bracken-{wildcards.feature_table_id}/meta.yaml"
    
    with open(meta_file, 'r') as f:
        meta = yaml.safe_load(f)
    
    # Reconstruct paths: kraken2/{qc_filter}/{sample}.bracken
    qc_filter = meta['qc_filter']
    samples = meta['samples']
    
    return [f"kraken2/{qc_filter}/{sample}.bracken" for sample in samples]


def get_abundance_metric(wildcards):
    """Read meta.yaml to get abundance metric."""
    meta_file = f"feature_tables/bracken-{wildcards.feature_table_id}/meta.yaml"
    
    with open(meta_file, 'r') as f:
        meta = yaml.safe_load(f)
    
    return meta.get('abundance_metric', 'new_est_reads')

rule build_bracken_abundance_table:
    input:
        files = get_bracken_table_inputs
    output:
        table = "feature_tables/bracken-{feature_table_id}/abundance_table.tsv"
    params:
        metric = get_abundance_metric,
        input_args = lambda wildcards, input: ' '.join([f'-i {f}' for f in input.files])
    shell:
        """
        python ~/MGX/scripts/build_abundance_table.py \
            {params.input_args} \
            --abundance-metric {params.metric} \
            --output {output.table}
        """

rule create_taxonomy_table:
    input:
        abundance = "feature_tables/bracken-{feature_table_id}/abundance_table.tsv"
    output:
        taxonomy = "feature_tables/bracken-{feature_table_id}/taxonomy_table.tsv"
    params:
        taxdump = DATABASES['taxdump']
    conda:
        "../../envs/taxonkit.yaml"
    shell:
        """
        # Extract tax IDs from first column (skip header, remove tax_id__ prefix)
        tail -n +2 {input.abundance} | cut -f1 | sed 's/tax_id__//' > {output.taxonomy}.tmp.ids
        
        # Get lineage and reformat with standard ranks
        cat {output.taxonomy}.tmp.ids | \
            taxonkit lineage --data-dir {params.taxdump} | \
            taxonkit reformat2 \
                --data-dir {params.taxdump} \
                -I 1 \
                -f "{{domain|acellular root|superkingdom}};{{phylum}};{{class}};{{order}};{{family}};{{genus}};{{species}}" \
                -r "NA" > {output.taxonomy}.tmp.reformatted
        
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
        """

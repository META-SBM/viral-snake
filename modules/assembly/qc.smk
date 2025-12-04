rule profile_contigs:
    """
    Generate per-contig statistics including length and GC content.
    Uses seqkit to extract contig-level metrics.
    """
    input:
        contigs = "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        stats = "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/contig_stats.tsv"
    wildcard_constraints:
        dataset = "[^/]+",
        prefix = "(assembly|co_assembly)/.+"
    conda:
        "../../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        echo -e "contig_id\tlength\tgc_content" > {output.stats}
        seqkit fx2tab --length --gc --name {input.contigs} >> {output.stats}
        """


rule summarize_contigs:
    """
    Generate assembly-level summary statistics.
    Provides overview metrics like N50, total length, contig count.
    """
    input:
        contigs = "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        summary = "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/contig_summary.tsv"
    wildcard_constraints:
        dataset = "[^/]+",
        prefix = "(assembly|co_assembly)/.+"
    conda:
        "../../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        seqkit stats -T {input.contigs} > {output.summary}
        """
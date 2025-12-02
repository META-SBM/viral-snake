rule profile_contigs:
    input:
        contigs = "{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        stats = "{prefix}/contigs_formatted_minlen_{min_len}/contig_stats.tsv"
    conda:
        "../../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        echo -e "contig_id\tlength\tgc_content" > {output.stats}
        seqkit fx2tab --length --gc --name {input.contigs} >> {output.stats}
        """
rule summarize_contigs:
    input:
        contigs = "{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        summary = "{prefix}/contigs_formatted_minlen_{min_len}/contig_summary.tsv"
    conda:
        "../../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        seqkit stats -T {input.contigs} > {output.summary}
        """
        

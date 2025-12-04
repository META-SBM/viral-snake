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
        
rule quast_individual:
    input:
        contigs = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        report = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/quast_report.html"
    params:
        outdir = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/",
        report = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/report.html"
    threads: 8
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        quast.py {input.contigs} \
            -o {params.outdir} \
            -t {threads} \
            --no-plots \
            --no-html \
            --no-icarus
        
        # Generate HTML report
        quast.py {input.contigs} \
            -o {params.outdir} \
            -t {threads} 

        mv {params.report} {output.report}
        """





rule quast_comparison:
    input:
        megahit = "assembly/megahit/{qc_filter}/{sample}/contigs.fa",
        metaspades = "assembly/metaspades/{qc_filter}/{sample}/contigs.fa"
    output:
        report = "qc/quast/comparison/{qc_filter}/{sample}/report.html"
    params:
        outdir = "qc/quast/comparison/{qc_filter}/{sample}"
    threads: 8
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        quast.py {input.megahit} {input.metaspades} \
            -o {params.outdir} \
            -t {threads} \
            --labels megahit,metaspades
        """
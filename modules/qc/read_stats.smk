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
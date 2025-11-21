rule subsample_reads:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        r1 = "reads/{qc_filter}__subsample_{n_reads}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}__subsample_{n_reads}/{sample}_R2.fastq.gz"
    params:
        seed = 42  # Random seed for reproducibility
    threads: 1
    log:
        "reads/{qc_filter}__subsample_{n_reads}/{sample}.log"
    conda:
        "../../envs/seqtk.yaml"
    shell:
        """
        seqtk sample -s {params.seed} {input.r1} {wildcards.n_reads} | gzip > {output.r1} 2> {log}
        seqtk sample -s {params.seed} {input.r2} {wildcards.n_reads} | gzip > {output.r2} 2>> {log}
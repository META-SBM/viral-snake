rule cutadapt:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        r1 = "reads/{qc_filter}__cutadapt_no_mgi_min_len_90/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}__cutadapt_no_mgi_min_len_90/{sample}_R2.fastq.gz"
    params:
        adapter_r1 = "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
        adapter_r2 = "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGT",
        quality = 10,
        min_length = 90
    threads: 8
    log:
        "reads/{qc_filter}__cutadapt_no_mgi_min_len_90/{sample}.log"
    conda:
        "../../envs/cutadapt.yaml"
    shell:
        """
        cutadapt \
            -a {params.adapter_r1} \
            -A {params.adapter_r2} \
            -q {params.quality} \
            -m {params.min_length} \
            --pair-filter=any \
            -j {threads} \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """


rule fastqc:
    input:
        "reads/{qc_filter}/{sample}_R{read}.fastq.gz"
    output:
        html = "qc/fastqc/{qc_filter}/{sample}_R{read}_fastqc.html",
        zip = "qc/fastqc/{qc_filter}/{sample}_R{read}_fastqc.zip"
    threads: 4
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o qc/fastqc/{wildcards.qc_filter}/ {input}
        """

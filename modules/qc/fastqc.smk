rule fastqc:
    """
    Run FastQC quality control on read files.
    Generates HTML report and zip archive with detailed metrics.
    """
    input:
        "{fs_prefix}/{dataset}/reads/{qc_filter}/{sample}_R{read}.fastq.gz"
    output:
        html = "{fs_prefix}/{dataset}/qc/fastqc/{qc_filter}/{sample}_R{read}_fastqc.html",
        zip = "{fs_prefix}/{dataset}/qc/fastqc/{qc_filter}/{sample}_R{read}_fastqc.zip"
    params:
        outdir = "{fs_prefix}/{dataset}/qc/fastqc/{qc_filter}"
    wildcard_constraints:
        dataset = "[^/]+",
        sample = "[^/]+",
        read = "[12]"
    threads: 4
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run FastQC
        fastqc -t {threads} -o {params.outdir} {input}
        """
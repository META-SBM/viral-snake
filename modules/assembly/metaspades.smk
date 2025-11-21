
rule metaspades:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        contigs = "assembly/metaspades/{qc_filter}/{sample}/contigs.fa"
    params:
        outdir = "assembly/metaspades/{qc_filter}/{sample}"
    threads: 32
    resources:
        mem_mb = 500000  # 500GB
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        rm -rf {params.outdir}
        metaspades.py -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            -t {threads} \
            -m {resources.mem_mb} \
            --only-assembler
        
        # metaSPAdes outputs contigs.fasta, rename to contigs.fa
        mv {params.outdir}/contigs.fasta {output.contigs}
        """

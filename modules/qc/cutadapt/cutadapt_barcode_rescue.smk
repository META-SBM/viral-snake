rule cutadapt_barcode_rescue_udi:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz",
        barcode_r1 = "config/barcodes_r1.fasta",
        barcode_r2 = "config/barcodes_r2.fasta"
    output:
        directory("reads/{qc_filter}__barcode_rescue/{sample}/")
    threads: 8
    log:
        "reads/{qc_filter}__barcode_rescue/{sample}.log"
    conda:
        "cutadapt.yaml"
    shell:
        """
        mkdir -p {output}
        
        cutadapt \
            -j {threads} \
            -e 0.15 \
            --no-indels \
            --pair-adapters \
            -a file:{input.barcode_r1} \
            -A file:{input.barcode_r2} \
            -o {output}/demultiplexed-{{name}}_R1.fastq.gz \
            -p {output}/demultiplexed-{{name}}_R2.fastq.gz \
            --discard-untrimmed \
            {input.r1} {input.r2} \
            > {log} 2>&1
        
        echo "Demultiplexing complete:" | tee -a {log}
        ls -lh {output}/ | tee -a {log}
        """
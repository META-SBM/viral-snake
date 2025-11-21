rule minimap2_align_to_references:
    """Minimap2 alignment of contigs to viral reference genomes"""
    input:
        contigs = "{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa",
        reference = "refseq_reference/{reference_org}/ncbi_dataset/data/genomic.fna"
    output:
        bam = "{prefix}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/alignments.sorted.bam",
        bai = "{prefix}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/alignments.sorted.bam.bai"
    params:
        preset = "asm5"  # Change to asm10 or asm20 for more divergent viruses
    threads: 6
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/minimap2.log"
    benchmark:
        "{prefix}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/minimap2.benchmark.txt"
    conda:
        "../../envs/minimap2.yaml"
    shell:
        """
        echo "=== Running Minimap2 alignment ===" | tee {log}
        echo "Input contigs: {input.contigs}" | tee -a {log}
        echo "Reference: {input.reference}" | tee -a {log}
        echo "Output BAM: {output.bam}" | tee -a {log}
        echo "Preset: {params.preset}" | tee -a {log}
        echo "Threads: {threads}" | tee -a {log}
        echo "" | tee -a {log}
        
        # Generate sorted BAM
        minimap2 -ax {params.preset} \
            -t {threads} \
            {input.reference} \
            {input.contigs} \
            2>> {log} | \
        samtools sort -@ {threads} -o {output.bam} - \
            2>> {log}
        
        # Index BAM
        samtools index -@ {threads} {output.bam} 2>> {log}
        
        echo "" | tee -a {log}
        echo "=== Minimap2 alignment complete ===" | tee -a {log}
        
        # Summary statistics
        num_bam_alignments=$(samtools view -c {output.bam})
        num_mapped=$(samtools view -F 4 -c {output.bam})
        echo "Total alignments: $num_bam_alignments" | tee -a {log}
        echo "Mapped contigs: $num_mapped" | tee -a {log}
        """

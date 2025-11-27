

def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(seq))

def get_sample_adapters(wildcards):
    """Generate R1 and R2 adapters for a specific sample"""
    # Get barcode number for this sample
    barcode_num = SAMPLE_TO_BARCODE_NUM[wildcards.sample]
    
    # Get 20bp index for this barcode
    index_20bp = BARCODE_NUM_TO_INDEX[barcode_num]
    
    # Generate R1 and R2 indexes
    index_r1 = index_20bp[-10:]  # Last 10bp for R1
    index_r2 = reverse_complement(index_20bp[:10])  # First 10bp revcomp for R2
    
    # Full adapters
    adapter_r1 = ADAPTER_R1_BASE + index_r1
    adapter_r2 = ADAPTER_R2_BASE + index_r2
    
    return {'r1': adapter_r1, 'r2': adapter_r2}


rule cutadapt_barcode_rescue:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        r1 = "reads/{qc_filter}__barcode_rescue/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}__barcode_rescue/{sample}_R2.fastq.gz",
        adapter_info = "reads/{qc_filter}__barcode_rescue/{sample}.adapters.yaml"
    params:
        adapters = get_sample_adapters,
        tmpdir = "reads/{qc_filter}__barcode_rescue/{sample}_tmp"
    threads: 4
    log:
        "reads/{qc_filter}__barcode_rescue/{sample}.log"
    conda:
        "cutadapt.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}

        # Save adapter configuration
        cat > {output.adapter_info} << EOF
sample: {wildcards.sample}
barcode_rescue_parameters:
  adapter_r1: "{params.adapters[r1]}"
  adapter_r2: "{params.adapters[r2]}"
  error_rate: 0.1
  no_indels: true
  concatenated_output: true
run_info:
  threads: {threads}
  timestamp: $(date -Iseconds)
EOF

        cutadapt \
            -e 0.1 \
            --no-indels \
            -a barcodeR1={params.adapters[r1]} \
            -A barcodeR2={params.adapters[r2]} \
            -o {params.tmpdir}/{{name1}}-{{name2}}_R1.fastq.gz \
            -p {params.tmpdir}/{{name1}}-{{name2}}_R2.fastq.gz \
            -j {threads} \
            {input.r1} {input.r2} \
            > {log} 2>&1
        
        echo "Sample: {wildcards.sample}" | tee -a {log}
        echo "R1 adapter: {params.adapters[r1]}" | tee -a {log}
        echo "R2 adapter: {params.adapters[r2]}" | tee -a {log}
        echo "" | tee -a {log}
        
        # Concatenate R1: barcoded + unknown
        cat {params.tmpdir}/barcodeR1-barcodeR2_R1.fastq.gz \
            {params.tmpdir}/unknown-unknown_R1.fastq.gz \
            > {output.r1}
        
        # Concatenate R2: barcoded + unknown
        cat {params.tmpdir}/barcodeR1-barcodeR2_R2.fastq.gz \
            {params.tmpdir}/unknown-unknown_R2.fastq.gz \
            > {output.r2}
        
        rm -rf {params.tmpdir}
        
        echo "Adapter configuration saved to: {output.adapter_info}" | tee -a {log}
        echo "Complete!" | tee -a {log}
        """
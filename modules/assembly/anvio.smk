def get_reformat_params(wildcards):
    """
    Extract proper prefix for contig naming with validation.
    """
    full_prefix = wildcards.prefix
    path_parts = full_prefix.split('/')
    
    # print(f"\n=== REFORMAT CONTIGS DEBUG ===")
    # print(f"Full prefix: {full_prefix}")
    
    # Parse based on expected structure
    if full_prefix.startswith("assembly/"):
        # Structure: assembly/{assembler}/{qc_filter}/{sample}
        if len(path_parts) >= 4:
            assembler = path_parts[1]
            qc_filter = path_parts[2]
            sample = path_parts[3]
            contig_prefix = sample
            # print(f"Type: Individual assembly")
            # print(f"Assembler: {assembler}")
            # print(f"QC filter: {qc_filter}")
            # print(f"Sample: {sample}")
        else:
            raise ValueError(f"Unexpected assembly path structure: {full_prefix}")
            
    elif full_prefix.startswith("co_assembly/"):
        # Structure: co_assembly/{assembler}/{collection}
        if len(path_parts) >= 3:
            assembler = path_parts[1]
            collection = path_parts[2]
            contig_prefix = collection
            # print(f"Type: Co-assembly")
            # print(f"Assembler: {assembler}")
            # print(f"Collection: {collection}")
        else:
            raise ValueError(f"Unexpected co_assembly path structure: {full_prefix}")
    else:
        raise ValueError(f"Unknown assembly type for prefix: {full_prefix}")
    
    # print(f"Contig prefix for anvio: {contig_prefix}")
    # print(f"Min length: {wildcards.min_len}")
    # print(f"=============================\n")
    
    return contig_prefix


rule reformat_contigs_unified:
    """Works for both individual assemblies and co-assemblies"""
    input:
        contigs = "{prefix}/contigs.fa"
    output:
        contigs = "{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa",
        report = "{prefix}/contigs_formatted_minlen_{min_len}/reformat_report.txt"
    params:
        contig_prefix = get_reformat_params,
        min_len = lambda wildcards: wildcards.min_len
    conda:
        "anvio-8"
    shell:
        """
        echo "=== Running anvi-script-reformat-fasta ==="
        echo "Input: {input.contigs}"
        echo "Output: {output.contigs}"
        echo "Prefix: {params.contig_prefix}"
        echo "Min length: {params.min_len}"
        echo ""
        
        anvi-script-reformat-fasta {input.contigs} \
            -o {output.contigs} \
            --simplify-names \
            --report-file {output.report} \
            --prefix {params.contig_prefix} \
            --min-len {params.min_len}
        
        echo ""
        echo "=== Reformatting complete ==="
        echo "Report saved to: {output.report}"
        """
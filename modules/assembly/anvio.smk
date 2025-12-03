def get_reformat_params(wildcards):
    """
    Extract proper prefix for contig naming with validation.
    
    Args:
        wildcards: Must contain 'prefix' and 'dataset'
        - prefix: assembly path like 'assembly/megahit/{qc}/{sample}' or 'co_assembly/megahit/{collection}'
        - dataset: dataset name for namespacing contigs
        
    Returns:
        str: Contig prefix in format '{dataset}_{sample}' or '{dataset}_{collection}'
    """
    full_prefix = wildcards.prefix
    dataset = wildcards.dataset
    path_parts = full_prefix.split('/')
    
    # Parse based on expected structure
    if full_prefix.startswith("assembly/"):
        # Structure: assembly/{assembler}/{qc_filter}/{sample}
        if len(path_parts) >= 4:
            assembler = path_parts[1]
            qc_filter = path_parts[2]
            sample = path_parts[3]
            # Include dataset in contig prefix to avoid collisions across datasets
            contig_prefix = f"{dataset}__{sample}"
        else:
            raise ValueError(f"Unexpected assembly path structure: {full_prefix}")
            
    elif full_prefix.startswith("co_assembly/"):
        # Structure: co_assembly/{assembler}/{collection}
        if len(path_parts) >= 3:
            assembler = path_parts[1]
            collection = path_parts[2]
            # Include dataset in contig prefix to avoid collisions across datasets
            contig_prefix = f"{dataset}__{collection}"
        else:
            raise ValueError(f"Unexpected co_assembly path structure: {full_prefix}")
    else:
        raise ValueError(f"Unknown assembly type for prefix: {full_prefix}")
    
    return contig_prefix


rule reformat_contigs_unified:
    """
    Reformat contig headers and filter by minimum length using anvi'o.
    Works for both individual assemblies and co-assemblies.
    Adds dataset prefix to contig names to prevent collisions across datasets.
    """
    input:
        contigs = "{fs_prefix}/{dataset}/{prefix}/contigs.fa"
    output:
        contigs = "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa",
        report = "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/reformat_report.txt"
    params:
        contig_prefix = get_reformat_params,
        min_len = lambda wildcards: wildcards.min_len
    wildcard_constraints:
        dataset = "[^/]+",
        prefix = "(assembly|co_assembly)/.+"
    conda:
        "anvio-8"
    shell:
        """
        echo "=== Running anvi-script-reformat-fasta ==="
        echo "Dataset: {wildcards.dataset}"
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
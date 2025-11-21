rule kraken2_classify_contigs:
    input:
        fasta = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        report = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/kraken2/confidence_{confidence}/kraken2.report",
        output = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/kraken2/confidence_{confidence}/kraken2.output"
    params:
        db = DATABASES['kraken2']
    log:
        "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/kraken2/confidence_{confidence}/kraken2.log"
    threads: THREADS['kraken2_contigs']
    conda:
        "../../envs/kraken2.yaml"
    shell:
        """
        k2 classify --db {params.db} \
            --threads {threads} \
            --confidence {wildcards.confidence} \
            --memory-mapping \
            --report {output.report} \
            --output {output.output} \
            --log {log} \
            {input.fasta}
        """

rule add_taxonomy_to_kraken2_contigs:
    input:
        kraken_output = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/kraken2/confidence_{confidence}/kraken2.output"
    output:
        kraken_taxonomy = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/kraken2/confidence_{confidence}/kraken2_output_with_taxonomy.tsv"
    params:
        taxdump = DATABASES['taxdump']
    conda:
        "../../envs/taxonkit.yaml"
    shell:
        """
        # Extract unique tax IDs from column 3 (skip unclassified with taxid 0)
        awk '$3 != "0" {{print $3}}' {input.kraken_output} | \
            sort -u > {output.kraken_taxonomy}.tmp.taxids
        
        # Get reformatted taxonomy using reformat2
        cat {output.kraken_taxonomy}.tmp.taxids | \
            taxonkit lineage --data-dir {params.taxdump} | \
            taxonkit reformat2 \
                --data-dir {params.taxdump} \
                -I 1 \
                -f "{{domain|acellular root|superkingdom}};{{phylum}};{{class}};{{order}};{{family}};{{genus}};{{species}}" \
                -r "NA" | \
            cut -f1,3 | \
            awk 'BEGIN {{FS="\t"; OFS="\t"}} 
                 {{
                   split($2, ranks, ";");
                   print $1, ranks[1], ranks[2], ranks[3], ranks[4], ranks[5], ranks[6], ranks[7]
                 }}' > {output.kraken_taxonomy}.tmp.lineages
        
        # Join kraken2 output with taxonomy
        awk 'BEGIN {{FS=OFS="\t"}}
             NR==FNR {{
               tax[$1] = $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8;
               next
             }}
             {{
               taxid = $3;
               if (taxid == "0") {{
                 # Unclassified sequences
                 print $1, $2, $3, $4, $5, "Unclassified\tNA\tNA\tNA\tNA\tNA\tNA"
               }} else if (taxid in tax) {{
                 print $1, $2, $3, $4, $5, tax[taxid]
               }} else {{
                 print $1, $2, $3, $4, $5, "NA\tNA\tNA\tNA\tNA\tNA\tNA"
               }}
             }}' {output.kraken_taxonomy}.tmp.lineages {input.kraken_output} > {output.kraken_taxonomy}.tmp.joined
        
        # Add header
        echo -e "Status\tSequenceID\tTaxID\tLength\tKmerInfo\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" | \
            cat - {output.kraken_taxonomy}.tmp.joined > {output.kraken_taxonomy}
        
        # Clean up
        rm {output.kraken_taxonomy}.tmp.*
        """


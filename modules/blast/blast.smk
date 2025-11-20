BLAST_DB = "/mnt/mgx/DATABASES/blast/core_nt/01-Nov-2025/core_nt"
rule blastn_nt:
    input:
        fasta = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        blast_out = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/blast/{database}/blastn.tsv"
    benchmark:
        "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/blast/{database}/blast.becnhmark.txt"
    params:
        db = BLAST_DB,
        max_target_seqs = 5,
        evalue = "1e-5",
        max_hsps = 3
    threads: 6
    conda:
        "blast"
    shell:
        """
        blastn -query {input.fasta} \
            -db {params.db} \
            -out {output.blast_out} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
            -num_threads {threads} \
            -max_target_seqs {params.max_target_seqs} \
            -max_hsps {params.max_hsps} \
            -evalue {params.evalue}
        """

rule add_taxonomy_to_blast:
    input:
        blast_out = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/blast/{database}/blastn.tsv"
    output:
        blast_taxonomy = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/blast/{database}/blastn_with_taxonomy.tsv"
    params:
        taxdump = "/mnt/mgx/DATABASES/taxdump/01-Nov-2025"
    conda:
        "taxonkit"
    shell:
        """
        # Extract unique tax IDs
        cut -f13 {input.blast_out} | \
            tr ';' '\n' | \
            grep -v '^$' | \
            grep -v 'N/A' | \
            sort -u > {output.blast_taxonomy}.tmp.taxids
        
        # Get reformatted taxonomy using reformat2 with pipe-separated rank alternatives
        cat {output.blast_taxonomy}.tmp.taxids | \
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
                 }}' > {output.blast_taxonomy}.tmp.lineages
        
        # Join BLAST results with taxonomy
        awk 'BEGIN {{FS=OFS="\t"}}
             NR==FNR {{
               tax[$1] = $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8;
               next
             }}
             {{
               taxid = $13;
               # Handle multiple taxids - take first one
               if (index(taxid, ";")) {{
                 split(taxid, ids, ";");
                 taxid = ids[1];
               }}
               if (taxid in tax) {{
                 print $0, tax[taxid]
               }} else {{
                 print $0, "NA\tNA\tNA\tNA\tNA\tNA\tNA"
               }}
             }}' {output.blast_taxonomy}.tmp.lineages {input.blast_out} > {output.blast_taxonomy}.tmp.joined
        
        # Add header
        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxids\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" | \
            cat - {output.blast_taxonomy}.tmp.joined > {output.blast_taxonomy}
        
        # Clean up
        rm {output.blast_taxonomy}.tmp.*
        """
rule calculate_contig_lca:
    """Calculate LCA for each contig using taxonkit"""
    input:
        hits = "{prefix}/hits_with_taxonomy.tsv"
    output:
        lca = "{prefix}/contig_lca.tsv",
        lca_detailed = "{prefix}/contig_lca_detailed.tsv"
    params:
        taxdump = "/mnt/mgx/DATABASES/taxdump/01-Nov-2025",
        top_n = 10
    log:
        "{prefix}/lca.log"
    conda:
        "taxonkit"
    shell:
        """
        echo "Starting LCA calculation" > {log}
        
        # Extract contig ID, bitscore, and taxid
        tail -n +2 {input.hits} | \
            awk -F'\t' 'BEGIN {{OFS="\t"}} {{print $1, $12, $14}}' | \
            sort -t$'\t' -k1,1 -k2,2rn > {output.lca}.tmp.sorted
        
        # Group by contig, take top N hits, collect taxids
        awk -F'\t' -v top_n={params.top_n} \
            'BEGIN {{OFS="\t"}}
             {{
                contig = $1;
                taxid = $3;
                
                # Handle multiple taxids - take first one
                if (index(taxid, ";")) {{
                    split(taxid, ids, ";");
                    taxid = ids[1];
                }}
                
                # Skip if no valid taxid
                if (taxid == "" || taxid == "NA") next;
                
                # Track count per contig
                if (contig != last_contig) {{
                    count = 0;
                    last_contig = contig;
                }}
                
                count++;
                if (count <= top_n) {{
                    if (contig in taxids) {{
                        taxids[contig] = taxids[contig] " " taxid;
                    }} else {{
                        taxids[contig] = taxid;
                    }}
                }}
             }}
             END {{
                for (contig in taxids) {{
                    print contig, taxids[contig];
                }}
             }}' {output.lca}.tmp.sorted > {output.lca}.tmp.taxids
        
        echo "Grouped into $(wc -l < {output.lca}.tmp.taxids) contigs" >> {log}
        
        # Calculate LCA using taxonkit
        cat {output.lca}.tmp.taxids | \
            taxonkit lca \
                --data-dir {params.taxdump} \
                -i 2 \
                -s " " \
                -D \
                -U \
                2>> {log} > {output.lca}.tmp.lca
        
        echo "Calculated LCA for $(wc -l < {output.lca}.tmp.lca) contigs" >> {log}
        
        # Extract just contig_id and lca_taxid
        cut -f1,3 {output.lca}.tmp.lca > {output.lca}.tmp.contig_taxid
        
        # Get taxid -> name mapping
        cut -f2 {output.lca}.tmp.contig_taxid | \
            taxonkit lineage --data-dir {params.taxdump} 2>> {log} | \
            awk -F'\t' '{{
                # Extract just the last element from lineage (the actual name)
                lineage = $2;
                n = split(lineage, parts, ";");
                name = parts[n];
                print $1, name;
            }}' > {output.lca}.tmp.names
        
        # Get formatted taxonomy
        cut -f2 {output.lca}.tmp.contig_taxid | \
            taxonkit reformat2 \
                --data-dir {params.taxdump} \
                -I 1 \
                -f "{{domain|acellular root|superkingdom}};{{phylum}};{{class}};{{order}};{{family}};{{genus}};{{species}}" \
                -r "NA" 2>> {log} | \
            cut -f2 > {output.lca}.tmp.formatted
        
        # Combine: contig_id, lca_taxid, lca_name, taxonomy_ranks
        paste {output.lca}.tmp.contig_taxid \
              <(cut -f2 {output.lca}.tmp.names) \
              {output.lca}.tmp.formatted | \
            awk 'BEGIN {{FS=OFS="\t"}}
                 {{
                   contig = $1;
                   lca_taxid = $2;
                   lca_name = $3;
                   taxonomy = $4;
                   
                   # Split taxonomy into ranks
                   split(taxonomy, ranks, ";");
                   
                   print contig, lca_taxid, lca_name, ranks[1], ranks[2], ranks[3], ranks[4], ranks[5], ranks[6], ranks[7];
                 }}' > {output.lca_detailed}.tmp
        
        # Add headers
        echo -e "contig_id\tlca_taxid\tlca_name\tlca_domain\tlca_phylum\tlca_class\tlca_order\tlca_family\tlca_genus\tlca_species" | \
            cat - {output.lca_detailed}.tmp > {output.lca_detailed}
        
        # Simple version
        echo -e "contig_id\tlca_taxid" > {output.lca}
        cut -f1,2 {output.lca_detailed}.tmp >> {output.lca}
        
        # Cleanup
        rm {output.lca}.tmp.*
        
        # Log summary
        echo "" >> {log}
        echo "=== LCA calculation complete ===" >> {log}
        echo "Total contigs: $(tail -n +2 {output.lca_detailed} | wc -l)" >> {log}
        echo "Top 5 LCA domains:" >> {log}
        tail -n +2 {output.lca_detailed} | cut -f4 | sort | uniq -c | sort -rn | head -5 >> {log}
        """
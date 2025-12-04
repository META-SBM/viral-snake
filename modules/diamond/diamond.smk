import yaml
from pathlib import Path

# Load DIAMOND presets
MODULE_DIR = Path(workflow.basedir) / "modules/diamond"
with open(MODULE_DIR / "diamond_presets.yaml") as f:
    DIAMOND_PRESETS = yaml.safe_load(f)['presets']

def get_diamond_params(wildcards):
    """Get all parameters for a DIAMOND preset"""
    preset = DIAMOND_PRESETS[wildcards.preset]
    return {
        'sensitivity_flag': preset['sensitivity_flag'],
        'block_size': preset['block_size'],
        'index_chunks': preset['index_chunks'],
        'max_target_seqs': preset['max_target_seqs'],
        'evalue': preset['evalue']
    }


rule diamond_blastx_unified:
    """DIAMOND BLASTX for both individual assemblies and co-assemblies"""
    input:
        contigs = "{prefix}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        hits = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/hits.txt"
    params:
        db = DATABASES['diamond_nr'],
        diamond_params = get_diamond_params
    threads: THREADS['diamond']
    wildcard_constraints:
        preset = "|".join(DIAMOND_PRESETS.keys())  # Validate preset names
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/diamond.log"
    benchmark:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/diamond.benchmark.txt"
    conda:
        "../../envs/diamond.yaml"
    wildcard_constraints:
        preset = "[^/]+"  # Explicit constraint: no slashes
    shell:
        """
        echo "=== Running DIAMOND BLASTX ===" | tee {log}
        echo "Input contigs: {input.contigs}" | tee -a {log}
        echo "Database: {params.db}" | tee -a {log}
        echo "Output: {output.hits}" | tee -a {log}
        echo "Preset: {wildcards.preset}" | tee -a {log}
        echo "Sensitivity: {params.diamond_params[sensitivity_flag]}" | tee -a {log}
        echo "Threads: {threads}" | tee -a {log}
        echo "" | tee -a {log}
        
        diamond blastx \
            -d {params.db} \
            -q {input.contigs} \
            -o {output.hits} \
            {params.diamond_params[sensitivity_flag]} \
            --block-size {params.diamond_params[block_size]} \
            --index-chunks {params.diamond_params[index_chunks]} \
            --threads {threads} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids \
            --max-target-seqs {params.diamond_params[max_target_seqs]} \
            --evalue {params.diamond_params[evalue]} \
            2>> {log}
        
        echo "" | tee -a {log}
        echo "=== DIAMOND BLASTX complete ===" | tee -a {log}
        num_hits=$(wc -l < {output.hits})
        echo "Total hits: $num_hits" | tee -a {log}
        """

rule add_taxonomy_to_diamond:
    input:
        diamond_hits = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/hits.txt"
    output:
        diamond_taxonomy = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/hits_with_taxonomy.tsv"
    params:
        taxdump = DATABASES['taxdump']
    conda:
        "../../envs/taxonkit.yaml"
    wildcard_constraints:
        preset = "[^/]+"  # Explicit constraint: no slashes
    shell:
        """
        # Extract unique tax IDs from column 14, handling multiple taxids separated by semicolons
        cut -f14 {input.diamond_hits} | \
            tr ';' '\n' | \
            grep -v '^$' | \
            sort -u > {output.diamond_taxonomy}.tmp.taxids
        
        # Get reformatted taxonomy using reformat2
        cat {output.diamond_taxonomy}.tmp.taxids | \
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
                 }}' > {output.diamond_taxonomy}.tmp.lineages
        
        # Join DIAMOND results with taxonomy
        # For hits with multiple taxids, use the first one for taxonomy lookup
        awk 'BEGIN {{FS=OFS="\t"}}
             NR==FNR {{
               tax[$1] = $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8;
               next
             }}
             {{
               taxid = $14;
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
             }}' {output.diamond_taxonomy}.tmp.lineages {input.diamond_hits} > {output.diamond_taxonomy}.tmp.joined
        
        # Add header
        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\tstaxids\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" | \
            cat - {output.diamond_taxonomy}.tmp.joined > {output.diamond_taxonomy}
        
        # Clean up
        rm {output.diamond_taxonomy}.tmp.*
        """
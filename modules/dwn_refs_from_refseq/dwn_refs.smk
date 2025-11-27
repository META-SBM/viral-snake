# ============================================================================
# modules/dwn_refs_from_refseq/dwn_refs.smk
# RefSeq genome downloading using checkpoints for dynamic targets
# ============================================================================

import pandas as pd
from pathlib import Path

MODULE_DIR = Path(workflow.basedir) / "modules/dwn_refs_from_refseq"


# ============================================================================
# Checkpoint: Extract taxid list from DIAMOND hits
# ============================================================================

checkpoint extract_species_for_download:
    """
    Extract unique taxids from filtered DIAMOND hits.
    This creates the list of genomes to download.
    If multiple taxids (semicolon-separated), takes the first one.
    """
    input:
        hits = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/hits.tsv"
    output:
        taxid_list = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/taxids_for_download.txt"
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/extract_taxids.log"
    shell:
        """
        # Extract taxid (col 14), split on semicolon and take first, get unique
        tail -n +2 {input.hits} | \
        awk -F'\t' '{{
            if ($14 != "") {{
                split($14, taxids, ";");
                print taxids[1]
            }}
        }}' | \
        sort -u \
        > {output.taxid_list} 2> {log}
        
        echo "Extracted $(wc -l < {output.taxid_list}) unique taxids" >> {log}
        """


# ============================================================================
# Input function: Read checkpoint output
# ============================================================================

def get_taxids_for_download(wildcards):
    """
    Read the checkpoint output and return list of taxids.
    Called after checkpoint execution completes.
    """
    # Access checkpoint output
    checkpoint_output = checkpoints.extract_species_for_download.get(**wildcards).output.taxid_list
    
    # Read the file
    taxids = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                taxids.append(line.strip())
    
    return taxids


# ============================================================================
# Dynamic rule: Download one genome per taxid
# ============================================================================

rule download_genome_by_taxid:
    """
    Download reference genome for a single taxonomy ID.
    One instance of this rule per taxid from checkpoint.
    """
    output:
        genome = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}/genomic.fna",
        metadata = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}/metadata.json"
    params:
        taxid = lambda w: w.taxid,
        outdir = lambda w, output: os.path.dirname(output.genome)
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}/download.log"
    threads: 1
    retries: 3
    resources:
        ncbi_downloads=1
    conda:
        "env.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Small random delay to stagger requests (0-5 seconds)
        sleep $((RANDOM % 5))
        
        # Download with timeout
        timeout 300 datasets download virus genome taxon {params.taxid} \
            --complete-only \
            --filename {params.outdir}/download.zip \
            2> {log} || {{
            echo "Download attempt failed for taxid {params.taxid}" >> {log}
        }}
        
        # Check if download succeeded
        if [ -f {params.outdir}/download.zip ] && [ -s {params.outdir}/download.zip ]; then
            # Extract
            unzip -q -o {params.outdir}/download.zip -d {params.outdir} 2>> {log}
            
            # Check if genome file exists
            if [ -f {params.outdir}/ncbi_dataset/data/genomic.fna ]; then
                # SUCCESS
                mv {params.outdir}/ncbi_dataset/data/genomic.fna {output.genome}
                
                genome_size=$(stat -f%z {output.genome} 2>/dev/null || stat -c%s {output.genome})
                seq_count=$(grep -c "^>" {output.genome})
                download_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
                
                cat > {output.metadata} <<EOF
{{
  "taxid": "{params.taxid}",
  "download_date": "$download_date",
  "genome_size": $genome_size,
  "sequence_count": $seq_count,
  "source": "NCBI RefSeq",
  "download_type": "complete",
  "status": "success"
}}
EOF
                
                echo "✓ Success: Downloaded taxid {params.taxid} ($seq_count sequences, $genome_size bytes)" >> {log}
            else
                # FAILED - no genome in archive
                echo "✗ Error: genomic.fna not found in archive for taxid {params.taxid}" >> {log}
                touch {output.genome}
                cat > {output.metadata} <<EOF
{{
  "taxid": "{params.taxid}",
  "download_date": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "status": "failed",
  "error": "genomic.fna not found in archive"
}}
EOF
            fi
        else
            # FAILED - download didn't work
            echo "✗ Error: Download failed for taxid {params.taxid}" >> {log}
            touch {output.genome}
            cat > {output.metadata} <<EOF
{{
  "taxid": "{params.taxid}",
  "download_date": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "status": "failed",
  "error": "Network error or no genomes available"
}}
EOF
        fi
        
        # Cleanup
        rm -rf {params.outdir}/ncbi_dataset {params.outdir}/download.zip {params.outdir}/README.md
        """


# ============================================================================
# Aggregate rule: Download all genomes
# ============================================================================

rule download_all_refseq_genomes:
    """
    Aggregate rule that triggers download of all species.
    Request this output to run the entire checkpoint workflow.
    """
    input:
        taxid_list = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/taxids_for_download.txt",
        genomes = lambda wildcards: expand(
            "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}/genomic.fna",
            prefix=wildcards.prefix,
            min_len=wildcards.min_len,
            preset=wildcards.preset,
            database=wildcards.database,
            filter_preset=wildcards.filter_preset,
            taxid=get_taxids_for_download(wildcards)
        ),
        metadata = lambda wildcards: expand(
            "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}/metadata.json",
            prefix=wildcards.prefix,
            min_len=wildcards.min_len,
            preset=wildcards.preset,
            database=wildcards.database,
            filter_preset=wildcards.filter_preset,
            taxid=get_taxids_for_download(wildcards)
        )
    output:
        summary = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/download_summary.txt"
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/download_summary.log"
    shell:
        """
        echo "RefSeq Download Summary" > {output.summary}
        echo "======================" >> {output.summary}
        echo "" >> {output.summary}
        
        # Count successes and failures
        success=0
        failed=0
        
        for metadata in {input.metadata}; do
            if grep -q '"status": "success"' $metadata 2>/dev/null; then
                ((success++))
            else
                ((failed++))
            fi
        done
        
        echo "Total taxids: $(echo {input.genomes} | wc -w)" >> {output.summary}
        echo "Successfully downloaded: $success" >> {output.summary}
        echo "Failed downloads: $failed" >> {output.summary}
        echo "" >> {output.summary}
        
        if [ $success -gt 0 ]; then
            echo "=== Successfully Downloaded ===" >> {output.summary}
            for metadata in {input.metadata}; do
                if grep -q '"status": "success"' $metadata 2>/dev/null; then
                    taxid=$(grep '"taxid"' $metadata | cut -d'"' -f4)
                    genome=$(dirname $metadata)/genomic.fna
                    size=$(du -h $genome 2>/dev/null | cut -f1)
                    echo "  ✓ Taxid $taxid ($size)" >> {output.summary}
                fi
            done
            echo "" >> {output.summary}
        fi
        
        if [ $failed -gt 0 ]; then
            echo "=== Failed Downloads ===" >> {output.summary}
            for metadata in {input.metadata}; do
                if ! grep -q '"status": "success"' $metadata 2>/dev/null; then
                    taxid=$(grep '"taxid"' $metadata | cut -d'"' -f4)
                    error=$(grep '"error"' $metadata | cut -d'"' -f4 || echo "Unknown error")
                    echo "  ✗ Taxid $taxid - $error" >> {output.summary}
                fi
            done
        fi
        
        cat {output.summary}
        """
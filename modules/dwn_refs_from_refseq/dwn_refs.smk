# ============================================================================
# modules/dwn_refs_from_refseq/dwn_refs.smk
# RefSeq genome downloading using checkpoints for dynamic targets
# ============================================================================

import pandas as pd
from pathlib import Path

MODULE_DIR = Path(workflow.basedir) / "modules/dwn_refs_from_refseq"


# ============================================================================
# Checkpoint: Extract taxid + species list from DIAMOND hits
# ============================================================================

checkpoint extract_species_for_download:
    """
    Extract unique taxids + species names from filtered DIAMOND hits.
    This creates the list of genomes to download.
    If multiple taxids (semicolon-separated), takes the first one.
    """
    input:
        hits = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/hits.tsv"
    output:
        taxid_list = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/taxids_for_download.tsv"
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/extract_taxids.log"
    run:
        import pandas as pd
        import re
        
        with open(log[0], 'w') as log_file:
            try:
                # Read DIAMOND hits
                df = pd.read_csv(input.hits, sep='\t')
                
                log_file.write(f"Read {len(df)} rows from DIAMOND hits\n")
                log_file.write(f"Columns: {list(df.columns)}\n\n")
                
                # Extract staxids and Species columns
                subset = df[['staxids', 'Species']].copy()
                
                # Drop rows with missing values
                subset = subset.dropna()
                subset = subset[subset['staxids'] != '']
                subset = subset[subset['Species'] != '']
                
                log_file.write(f"After filtering empty values: {len(subset)} rows\n")
                
                # Handle semicolon-separated taxids (take first)
                subset['staxids'] = subset['staxids'].astype(str).str.split(';').str[0]
                
                # Sanitize species names
                def sanitize_species(name):
                    """Remove problematic characters and sanitize for filesystem"""
                    name = str(name)
                    # Remove: ' " ` \ / ( ) [ ] { } and other problematic chars
                    name = re.sub(r"['\"`.:/\\()\[\]{}]", '', name)
                    # Replace whitespace with underscore
                    name = re.sub(r'\s+', '_', name)
                    # Collapse multiple underscores
                    name = re.sub(r'_+', '_', name)
                    # Remove leading/trailing underscores
                    name = name.strip('_')
                    return name
                
                subset['Species'] = subset['Species'].apply(sanitize_species)
                
                # Remove duplicates
                subset = subset.drop_duplicates()
                
                log_file.write(f"After removing duplicates: {len(subset)} unique taxid-species pairs\n\n")
                
                # Rename columns for output
                subset.columns = ['taxid', 'species']
                
                # Write to output
                subset.to_csv(output.taxid_list, sep='\t', index=False)
                
                log_file.write(f"✓ Extracted {len(subset)} unique taxids\n")
                log_file.write("✓ Species names sanitized for filesystem safety\n\n")
                log_file.write(f"First 10 entries:\n")
                log_file.write(subset.head(10).to_string(index=False))
                log_file.write("\n")
                
            except Exception as e:
                log_file.write(f"ERROR: {str(e)}\n")
                import traceback
                log_file.write(traceback.format_exc())
                raise


# ============================================================================
# Input function: Read checkpoint output
# ============================================================================

def get_taxids_and_species_for_download(wildcards):
    """
    Read the checkpoint output and return dict with taxids and species.
    Called after checkpoint execution completes.
    """
    # Access checkpoint output
    checkpoint_output = checkpoints.extract_species_for_download.get(**wildcards).output.taxid_list
    
    # Read the TSV file
    taxids = []
    species = []
    
    with open(checkpoint_output) as f:
        next(f)  # Skip header
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    taxid, sp = parts
                    # Replace spaces with underscores for folder names
                    species_safe = sp.replace(' ', '_')
                    taxids.append(taxid)
                    species.append(species_safe)
    
    return {'taxids': taxids, 'species': species}


# ============================================================================
# Dynamic rule: Download one genome per taxid
# ============================================================================

rule download_genome_by_taxid:
    """
    Download reference genome for a single taxonomy ID.
    One instance of this rule per taxid from checkpoint.
    Folder name includes both taxid and species for clarity.
    """
    input:
        taxid_list = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/taxids_for_download.tsv"
    output:
        genome = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}__{species}/genomic.fna",
        protein = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}__{species}/protein.faa",
        metadata = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}__{species}/metadata.json",
        archive = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}__{species}/download.zip"
    params:
        taxid = lambda w: w.taxid,
        species = lambda w: w.species.replace('_', ' '),
        species_safe = lambda w: w.species,
        outdir = lambda w, output: os.path.dirname(output.genome)
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}__{species}/download.log"
    threads: 1
    retries: 3
    resources:
        ncbi_downloads=1
    conda:
        "env.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        echo "Downloading taxid {params.taxid} ({params.species})" > {log}
        
        # Small random delay to stagger requests (0-5 seconds)
        sleep $((RANDOM % 5))
        
        # Download with timeout - include everything
        timeout 300 datasets download virus genome taxon {params.taxid} \
            --include genome,protein,cds,annotation \
            --filename {output.archive} \
            2>> {log} || {{
            echo "Download attempt failed for taxid {params.taxid} ({params.species})" >> {log}
        }}
        
        # Check if download succeeded
        if [ -f {output.archive} ] && [ -s {output.archive} ]; then
            # Extract
            unzip -q -o {output.archive} -d {params.outdir} 2>> {log}
            
            # Check if data directory exists
            if [ -d {params.outdir}/ncbi_dataset/data ]; then
                # Move ALL files from ncbi_dataset/data/ to output directory
                echo "Moving all files from ncbi_dataset/data/..." >> {log}
                mv {params.outdir}/ncbi_dataset/data/* {params.outdir}/ 2>> {log}
                
                # List what we got
                echo "Files extracted:" >> {log}
                ls -lh {params.outdir}/*.* >> {log} 2>&1
                
                # Check if required genome file exists
                if [ -f {output.genome} ]; then
                    # SUCCESS - create metadata
                    genome_size=$(stat -f%z {output.genome} 2>/dev/null || stat -c%s {output.genome})
                    seq_count=$(grep -c "^>" {output.genome})
                    
                    # Check what files we have
                    has_protein="false"
                    has_cds="false"
                    has_annotation="false"
                    protein_count=0
                    
                    if [ -s {output.protein} ]; then
                        has_protein="true"
                        protein_count=$(grep -c "^>" {output.protein})
                    else
                        touch {output.protein}
                    fi
                    
                    if [ -f {params.outdir}/cds.fna ] && [ -s {params.outdir}/cds.fna ]; then
                        has_cds="true"
                    fi
                    
                    if [ -f {params.outdir}/annotation_report.jsonl ]; then
                        has_annotation="true"
                    fi
                    
                    download_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
                    
                    cat > {output.metadata} <<EOF
{{
  "taxid": "{params.taxid}",
  "species": "{params.species}",
  "download_date": "$download_date",
  "genome_size": $genome_size,
  "sequence_count": $seq_count,
  "has_protein": $has_protein,
  "protein_count": $protein_count,
  "has_cds": $has_cds,
  "has_annotation": $has_annotation,
  "source": "NCBI RefSeq",
  "download_type": "complete",
  "status": "success"
}}
EOF
                    
                    echo "✓ Success: Downloaded taxid {params.taxid} ({params.species})" >> {log}
                    echo "  - Genome sequences: $seq_count" >> {log}
                    echo "  - Protein sequences: $protein_count" >> {log}
                    echo "  - Size: $genome_size bytes" >> {log}
                else
                    # FAILED - no genome file
                    echo "✗ Error: genomic.fna not found for taxid {params.taxid} ({params.species})" >> {log}
                    touch {output.genome}
                    touch {output.protein}
                    cat > {output.metadata} <<EOF
{{
  "taxid": "{params.taxid}",
  "species": "{params.species}",
  "download_date": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "status": "failed",
  "error": "genomic.fna not found in archive"
}}
EOF
                fi
            else
                # FAILED - no data directory
                echo "✗ Error: ncbi_dataset/data directory not found" >> {log}
                touch {output.genome}
                touch {output.protein}
                cat > {output.metadata} <<EOF
{{
  "taxid": "{params.taxid}",
  "species": "{params.species}",
  "download_date": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "status": "failed",
  "error": "ncbi_dataset/data directory not found in archive"
}}
EOF
            fi
        else
            # FAILED - download didn't work
            echo "✗ Error: Download failed for taxid {params.taxid} ({params.species})" >> {log}
            touch {output.genome}
            touch {output.protein}
            touch {output.archive}
            cat > {output.metadata} <<EOF
{{
  "taxid": "{params.taxid}",
  "species": "{params.species}",
  "download_date": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "status": "failed",
  "error": "Network error or no genomes available"
}}
EOF
        fi
        
        # Cleanup only the empty ncbi_dataset folder structure
        rm -rf {params.outdir}/ncbi_dataset {params.outdir}/README.md
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
        taxid_list = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/taxids_for_download.tsv",
        genomes = lambda wildcards: expand(
            "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}__{species}/genomic.fna",
            zip,
            prefix=[wildcards.prefix] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            min_len=[wildcards.min_len] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            preset=[wildcards.preset] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            database=[wildcards.database] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            filter_preset=[wildcards.filter_preset] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            taxid=get_taxids_and_species_for_download(wildcards)['taxids'],
            species=get_taxids_and_species_for_download(wildcards)['species']
        ),
        metadata = lambda wildcards: expand(
            "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/taxid_{taxid}__{species}/metadata.json",
            zip,
            prefix=[wildcards.prefix] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            min_len=[wildcards.min_len] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            preset=[wildcards.preset] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            database=[wildcards.database] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            filter_preset=[wildcards.filter_preset] * len(get_taxids_and_species_for_download(wildcards)['taxids']),
            taxid=get_taxids_and_species_for_download(wildcards)['taxids'],
            species=get_taxids_and_species_for_download(wildcards)['species']
        )
    output:
        summary = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/download_summary.txt"
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/refseq/download_summary.log"
    shell:
        """
        set +e  # Don't exit on errors in this script
        
        echo "RefSeq Download Summary" > {output.summary}
        echo "======================" >> {output.summary}
        echo "" >> {output.summary}
        
        # Count successes and failures
        success=0
        failed=0
        
        for metadata in {input.metadata}; do
            if [ -f "$metadata" ]; then
                if grep -q '"status": "success"' "$metadata" 2>/dev/null; then
                    success=$((success + 1))
                else
                    failed=$((failed + 1))
                fi
            fi
        done
        
        total=$(echo {input.genomes} | wc -w)
        
        echo "Total taxids: $total" >> {output.summary}
        echo "Successfully downloaded: $success" >> {output.summary}
        echo "Failed downloads: $failed" >> {output.summary}
        echo "" >> {output.summary}
        
        if [ $success -gt 0 ]; then
            echo "=== Successfully Downloaded ===" >> {output.summary}
            for metadata in {input.metadata}; do
                if [ -f "$metadata" ] && grep -q '"status": "success"' "$metadata" 2>/dev/null; then
                    taxid=$(grep '"taxid"' "$metadata" | cut -d'"' -f4)
                    species=$(grep '"species"' "$metadata" | cut -d'"' -f4)
                    genome=$(dirname "$metadata")/genomic.fna
                    if [ -f "$genome" ]; then
                        size=$(du -h "$genome" 2>/dev/null | cut -f1)
                        echo "  ✓ [$taxid] $species ($size)" >> {output.summary}
                    fi
                fi
            done
            echo "" >> {output.summary}
        fi
        
        if [ $failed -gt 0 ]; then
            echo "=== Failed Downloads ===" >> {output.summary}
            for metadata in {input.metadata}; do
                if [ -f "$metadata" ]; then
                    if ! grep -q '"status": "success"' "$metadata" 2>/dev/null; then
                        taxid=$(grep '"taxid"' "$metadata" | cut -d'"' -f4)
                        species=$(grep '"species"' "$metadata" | cut -d'"' -f4 2>/dev/null || echo "Unknown")
                        error=$(grep '"error"' "$metadata" | cut -d'"' -f4 2>/dev/null || echo "Unknown error")
                        echo "  ✗ [$taxid] $species - $error" >> {output.summary}
                    fi
                fi
            done
        fi
        
        echo "" >> {log}
        echo "Summary created successfully" >> {log}
        cat {output.summary} >> {log}
        cat {output.summary}
        """
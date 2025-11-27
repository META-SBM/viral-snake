import yaml
from pathlib import Path

# Load diamond filter presets
MODULE_DIR = Path(workflow.basedir) / "modules/diamond/filter"

with open(MODULE_DIR / "presets.yaml") as f:
    DIAMOND_FILTER_PRESETS = yaml.safe_load(f)['filters']


def build_filter_command(wildcards):
    """Build filter command from preset configuration"""
    preset = DIAMOND_FILTER_PRESETS[wildcards.filter_preset]
    
    cmd_parts = [
        f"--min-pident {preset.get('min_pident', 50)}",
        f"--min-length {preset.get('min_length', 100)}",
    ]
    
    if 'min_bitscore' in preset:
        cmd_parts.append(f"--min-bitscore {preset['min_bitscore']}")
    
    for domain in preset.get('keep_domains', []):
        cmd_parts.append(f"--keep-domains '{domain}'")
    
    for domain in preset.get('exclude_domains', []):
        cmd_parts.append(f"--exclude-domains '{domain}'")
    
    return " ".join(cmd_parts)


rule filter_diamond:
    """Filter DIAMOND results using predefined presets"""
    input:
        hits = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/hits_with_taxonomy.tsv"
    output:
        hits = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/hits.tsv",
        summary = "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/summary.txt"
    params:
        filter_cmd = build_filter_command
    wildcard_constraints:
        filter_preset = "|".join(DIAMOND_FILTER_PRESETS.keys())
    log:
        "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/{filter_preset}/filter.log"
    shell:
        """
        python {MODULE_DIR}/filter_diamond.py \
            {input.hits} \
            {output.hits} \
            {params.filter_cmd} \
            2> {log}
        
        # Create summary
        echo "Filter preset: {wildcards.filter_preset}" > {output.summary}
        echo "Parameters: {params.filter_cmd}" >> {output.summary}
        echo "Input hits: $(wc -l < {input.hits})" >> {output.summary}
        echo "Output hits: $(wc -l < {output.hits})" >> {output.summary}
        """



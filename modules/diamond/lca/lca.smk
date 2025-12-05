import yaml
from pathlib import Path

# Load diamond filter presets
LCA_MODULE_DIR = Path(workflow.basedir) / "modules/diamond/lca"

rule lca:
    input: 
        hits = "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/hits_with_taxonomy.tsv"

    output: 
        "{fs_prefix}/{dataset}/{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/LCA.tsv"

    params: 
        coef = 0.99,
        taxdump = DATABASES['taxdump']

    conda:
        "../../../envs/taxonkit.yaml"
    wildcard_constraints:
        preset = "[^/]+",
        dataset = "[^/]+",
        prefix = "(assembly|co_assembly)/.+"
    shell:
         "python3 {LCA_MODULE_DIR}/lca.py -input {input.hits} -output {output} -coef {params.coef} -taxdump {params.taxdump}"
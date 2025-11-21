import yaml
from pathlib import Path

# Load diamond filter presets
MODULE_DIR = Path(workflow.basedir) / "modules/lca"

rule lca:
    input: "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/hits_with_taxonomy.tsv"

    output: "{prefix}/contigs_formatted_minlen_{min_len}/diamond_{preset}/{database}/LCA.tsv"

    params: 
        coef = 0.99

    conda:
        "../../envs/taxonkit.yaml"
    wildcard_constraints:
        preset = "[^/]+"  # Explicit constraint: no slashes

    shell:
         "python3 {MODULE_DIR}/lca.py -input {input} -output {output} -coef {params.coef}"


from pathlib import Path

# Path to plotting script
MOSDEPTH_PLOT_DIR = Path(workflow.basedir) / "modules/mosdepth"
PLOT_SCRIPT = MOSDEPTH_PLOT_DIR / "plot_mosdepth.py"


rule plot_mosdepth:
    """
    Create interactive coverage plots from mosdepth output.
    Generates both HTML and pickle formats.
    """
    input:
        regions = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/depth.mosdepth__window_{window_size}.regions.bed.gz"
    output:
        html = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/depth.mosdepth__window_{window_size}.plot.html",
        pickle = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/depth.mosdepth__window_{window_size}.plot.pickle"
    params:
        script = PLOT_SCRIPT,
        sample_name = lambda wildcards: wildcards.sample
    wildcard_constraints:
        dataset = "[^/]+",
        reference_path = ".+",
        query_qc = "[^/]+",
        sample = "[^/]+",
        window_size = "\d+"
    log:
        "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/depth.mosdepth__window_{window_size}.plot.log"
    conda:
        "env.yaml"
    shell:
        """
        python {params.script} \
            {input.regions} \
            {output.html} \
            {output.pickle} \
            --sample-name {params.sample_name} \
            2> {log}
        """
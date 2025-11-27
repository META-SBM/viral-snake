import yaml
from pathlib import Path

# Load cutadapt presets
MODULE_DIR = Path(workflow.basedir) / "modules/qc/cutadapt"
with open(MODULE_DIR / "cutadapt_presets.yaml") as f:
    CUTADAPT_PRESETS = yaml.safe_load(f)['presets']

def get_cutadapt_params(wildcards):
    """Get all parameters for a cutadapt preset"""
    preset = CUTADAPT_PRESETS[wildcards.cutadapt_preset]
    return {
        'adapter_r1': preset['adapter_r1'],
        'adapter_r2': preset['adapter_r2'],
        'quality': preset['quality'],
        'min_length': preset['min_length'],
        'discard_trimmed': preset['discard_trimmed'],
        'pair_filter': preset['pair_filter'],
        'description': preset['description']
    }

rule cutadapt:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        r1 = "reads/{qc_filter}__cutadapt_{cutadapt_preset}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}__cutadapt_{cutadapt_preset}/{sample}_R2.fastq.gz",
        preset_record = "reads/{qc_filter}__cutadapt_{cutadapt_preset}/{sample}.preset.yaml"
    params:
        cutadapt_params = get_cutadapt_params
    threads: THREADS['cutadapt']
    wildcard_constraints:
        cutadapt_preset = "|".join(CUTADAPT_PRESETS.keys())  # Validate preset names
    log:
        "reads/{qc_filter}__cutadapt_{cutadapt_preset}/{sample}.log"
    benchmark:
        "reads/{qc_filter}__cutadapt_{cutadapt_preset}/{sample}.benchmark.txt"
    conda:
        "cutadapt.yaml"
    shell:
        """
        # Save preset configuration for reproducibility
        cat > {output.preset_record} << 'EOF'
preset_name: {wildcards.cutadapt_preset}
preset_description: "{params.cutadapt_params[description]}"
parameters:
  adapter_r1: "{params.cutadapt_params[adapter_r1]}"
  adapter_r2: "{params.cutadapt_params[adapter_r2]}"
  quality: {params.cutadapt_params[quality]}
  min_length: {params.cutadapt_params[min_length]}
  discard_trimmed: {params.cutadapt_params[discard_trimmed]}
  pair_filter: "{params.cutadapt_params[pair_filter]}"
run_info:
  sample: {wildcards.sample}
  threads: {threads}
  timestamp: $(date -Iseconds)
EOF
        
        echo "=== Running Cutadapt ===" | tee {log}
        echo "Input R1: {input.r1}" | tee -a {log}
        echo "Input R2: {input.r2}" | tee -a {log}
        echo "Output R1: {output.r1}" | tee -a {log}
        echo "Output R2: {output.r2}" | tee -a {log}
        echo "Preset: {wildcards.cutadapt_preset}" | tee -a {log}
        echo "Quality threshold: {params.cutadapt_params[quality]}" | tee -a {log}
        echo "Minimum length: {params.cutadapt_params[min_length]}" | tee -a {log}
        echo "Threads: {threads}" | tee -a {log}
        echo "" | tee -a {log}
        
        # Build cutadapt command with conditional parameters
        ADAPTER_R1="{params.cutadapt_params[adapter_r1]}"
        ADAPTER_R2="{params.cutadapt_params[adapter_r2]}"
        DISCARD_TRIMMED="{params.cutadapt_params[discard_trimmed]}"
        
        ADAPTER_ARGS=""
        if [ -n "$ADAPTER_R1" ]; then
            ADAPTER_ARGS="$ADAPTER_ARGS -a $ADAPTER_R1"
        fi
        if [ -n "$ADAPTER_R2" ]; then
            ADAPTER_ARGS="$ADAPTER_ARGS -A $ADAPTER_R2"
        fi
        
        DISCARD_ARGS=""
        if [ "$DISCARD_TRIMMED" = "True" ]; then
            DISCARD_ARGS="--discard-trimmed"
        fi
        
        cutadapt \
            $ADAPTER_ARGS \
            $DISCARD_ARGS \
            -q {params.cutadapt_params[quality]} \
            -m {params.cutadapt_params[min_length]} \
            --pair-filter={params.cutadapt_params[pair_filter]} \
            -j {threads} \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2} \
            >> {log} 2>&1
        
        echo "" | tee -a {log}
        echo "=== Cutadapt complete ===" | tee -a {log}
        echo "Preset configuration saved to: {output.preset_record}" | tee -a {log}
        """
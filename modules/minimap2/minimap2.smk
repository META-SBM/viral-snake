# viral-snake/modules/minimap2/minimap2.smk
# Minimap2 for assembly-to-reference alignment

# ============================================================================
# Helper Functions (copied from strobealign pattern)
# ============================================================================

def find_reference_file(ref_base):
    """
    Find reference file with any common FASTA extension.
    """
    extensions = ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz']
    for ext in extensions:
        fasta = f"{ref_base}{ext}"
        if os.path.exists(fasta):
            return fasta, ext
    return None, None


def get_reference_fasta_minimap2(wildcards):
    """
    Intelligent reference selection for minimap2.
    Supports: references/, assembly/, co_assembly/, refseq/
    """
    ref_path = wildcards.reference_path
    
    # User-provided references
    if ref_path.startswith("references/"):
        if ".sanitized" in ref_path:
            ref_base = ref_path.replace(".sanitized", "")
            full_base = f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_base}"
            original, ext = find_reference_file(full_base)
            if not original:
                raise FileNotFoundError(f"No FASTA file found for {full_base}")
            return f"{full_base}.sanitized{ext}"
        else:
            full_base = f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}"
            original, ext = find_reference_file(full_base)
            if not original:
                raise FileNotFoundError(f"No FASTA file found for {full_base}")
            return original
    
    # Use another assembly as reference!
    elif ref_path.startswith(("assembly/", "co_assembly/")):
        return f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}/contigs.fa"
    
    # RefSeq downloads
    elif ref_path.startswith("refseq/"):
        return f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}/genomic.fna"
    
    else:
        # Default: try as assembly path
        return f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}/contigs.fa"


# ============================================================================
# Rules
# ============================================================================

rule minimap2_align_assembly:
    """
    Align assembled contigs to a reference genome using minimap2.
    Uses asm5 preset for ~5% divergence (assembly vs reference).
    
    Reference can be:
    - references/my_genome (user reference, any extension)
    - refseq/GCF_000001405 (RefSeq download)
    - assembly/sample1/spades (another assembly!)
    - co_assembly/all/megahit (co-assembly as reference)
    """
    input:
        ref = get_reference_fasta_minimap2,  # INTELLIGENT!
        contigs = "{fs_prefix}/{dataset}/{query_path}/contigs.fa"
    output:
        bam = "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/alignments.sorted.bam",
        bai = "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/alignments.sorted.bam.bai"
    params:
        preset = "asm5"
    threads: THREADS.get('minimap2', 16)
    wildcard_constraints:
        dataset = "[^/]+",
        reference_path = ".+",
        query_path = r"(assembly|co_assembly)/.+"
    log:
        "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/alignment.log"
    benchmark:
        "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/benchmark.txt"
    conda:
        "../../envs/minimap2.yaml"
    shell:
        """
        echo "=== Minimap2 assembly alignment ===" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        echo "Query: {wildcards.query_path}" >> {log}
        echo "Reference path: {wildcards.reference_path}" >> {log}
        echo "Reference file: {input.ref}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "" >> {log}
        
        mkdir -p $(dirname {output.bam})
        
        minimap2 \
            -ax {params.preset} \
            -t {threads} \
            {input.ref} \
            {input.contigs} \
            2>> {log} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam} \
            - \
            2>> {log}
        
        samtools index -@ {threads} {output.bam} 2>> {log}
        
        echo "" >> {log}
        echo "=== Alignment complete ===" >> {log}
        """


rule minimap2_assembly_coverage:
    """
    Calculate how much of the reference is covered by assembly contigs.
    """
    input:
        bam = "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/alignments.sorted.bam",
        bai = "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/alignments.sorted.bam.bai"
    output:
        coverage = "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/reference_coverage.tsv"
    threads: 4
    wildcard_constraints:
        dataset = "[^/]+",
        reference_path = ".+",
        query_path = r"(assembly|co_assembly)/.+"
    log:
        "{fs_prefix}/{dataset}/alignment/minimap2/{reference_path}/__contigs__/{query_path}/coverage.log"
    conda:
        "../../envs/minimap2.yaml"
    shell:
        """
        echo "Calculating reference coverage" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        
        samtools coverage {input.bam} > {output.coverage} 2>> {log}
        
        echo "Coverage calculation complete!" >> {log}
        """
# viral-snake/modules/strobealign/strobealign.smk
# Strobealign alignment module
# Maps reads to assemblies, co-assemblies, or reference genomes

# ============================================================================
# Configuration
# ============================================================================

# Default read length for index creation (Illumina standard)
STROBEALIGN_READ_LENGTH = config.get('strobealign', {}).get('canonical_read_length', 150)

# ============================================================================
# Helper Functions
# ============================================================================

def find_reference_file(ref_base):
    """
    Find reference file with any common FASTA extension.
    
    Args:
        ref_base: Base path without extension
        
    Returns:
        tuple: (full_path, extension) or (None, None) if not found
    """
    extensions = ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz']
    for ext in extensions:
        fasta = f"{ref_base}{ext}"
        if os.path.exists(fasta):
            return fasta, ext
    return None, None


def get_reference_fasta(wildcards):
    """
    Construct path to reference FASTA from wildcards.
    
    Args:
        wildcards: Must contain fs_prefix, dataset, and reference_path
        
    Returns:
        str: Full path to reference FASTA file
    """
    ref_path = wildcards.reference_path
    
    # For references directory, handle extension detection
    if ref_path.startswith("references/"):
        # Check if sanitized version is requested
        if ".sanitized" in ref_path:
            # Extract base path (remove .sanitized)
            ref_base = ref_path.replace(".sanitized", "")
            full_base = f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_base}"
            original, ext = find_reference_file(full_base)
            if not original:
                raise FileNotFoundError(f"No FASTA file found for {full_base}")
            return f"{full_base}.sanitized{ext}"
        else:
            # Return original file
            full_base = f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}"
            original, ext = find_reference_file(full_base)
            if not original:
                raise FileNotFoundError(f"No FASTA file found for {full_base}")
            return original
    
    # For assemblies/co-assemblies
    elif ref_path.startswith(("assembly/", "co_assembly/")):
        fasta = f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}/contigs.fa"
    
    # For refseq downloads
    elif ref_path.startswith("refseq/"):
        fasta = f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}/genomic.fna"
    
    else:
        # Default: assume it's an assembly path
        fasta = f"{wildcards.fs_prefix}/{wildcards.dataset}/{ref_path}/contigs.fa"
    
    # Validation
    # if not os.path.exists(fasta):
    #     raise FileNotFoundError(f"Reference FASTA not found: {fasta}")
    
    return fasta


rule sanitize_reference_headers:
    """
    Replace spaces in FASTA headers with ___.
    Works with any FASTA extension (.fa, .fasta, .fna, .gz variants).
    """
    input:
        ref = "{fs_prefix}/{dataset}/references/{ref_path}.{ext}"
    output:
        sanitized = "{fs_prefix}/{dataset}/references/{ref_path}.sanitized.{ext}"
    wildcard_constraints:
        dataset = "[^/]+",
        ext = r"(fasta|fa|fna)(\.gz)?"
    log:
        "{fs_prefix}/{dataset}/references/{ref_path}.sanitized.{ext}.log"
    conda:
        "env.yaml"
    shell:
        """
        echo "Sanitizing headers in {input.ref}" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        
        # Handle compressed and uncompressed files
        if [[ {input.ref} == *.gz ]]; then
            zcat {input.ref} | sed '/^>/s/ /___/g' | gzip > {output.sanitized}
        else
            sed '/^>/s/ /___/g' {input.ref} > {output.sanitized}
        fi
        
        echo "Created: {output.sanitized}" >> {log}
        """


def get_strobealign_index(wildcards):
    """
    Get path to strobealign index file.
    
    Args:
        wildcards: Must contain fs_prefix, dataset, and reference_path
        
    Returns:
        str: Path to index file
    """
    return f"{wildcards.fs_prefix}/{wildcards.dataset}/alignment/__index__/strobealign/{wildcards.reference_path}/r{STROBEALIGN_READ_LENGTH}.sti"


rule strobealign_index:
    """
    Create strobealign index for a reference.
    The index is read-length optimized.
    """
    input:
        ref = get_reference_fasta
    output:
        index = "{fs_prefix}/{dataset}/alignment/__index__/strobealign/{reference_path}/r{read_len}.sti"
    params:
        read_len = lambda wildcards: wildcards.read_len
    threads: 1
    wildcard_constraints:
        dataset = "[^/]+",
        reference_path = ".+"
    log:
        "{fs_prefix}/{dataset}/alignment/__index__/strobealign/{reference_path}/r{read_len}.log"
    conda:
        "env.yaml"
    shell:
        """
        echo "Creating strobealign index" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        echo "Reference path: {wildcards.reference_path}" >> {log}
        echo "Read length: {params.read_len}" >> {log}
        echo "Reference file: {input.ref}" >> {log}
        
        strobealign --create-index \
            -r {params.read_len} \
            "{input.ref}" \
            2>> {log}
        
        # Create output directory
        mkdir -p $(dirname {output.index})
        
        # Move index to organized location
        # Strobealign creates index as: reference.fa.rXXX.sti
        mv "{input.ref}.r{params.read_len}.sti" {output.index}
        
        echo "Index created: {output.index}" >> {log}
        """


rule strobealign_align:
    """
    Align reads to reference using strobealign.
    Produces sorted BAM with index.
    """
    input:
        ref = get_reference_fasta,
        index = get_strobealign_index,
        r1 = "{fs_prefix}/{dataset}/reads/{query_qc}/{sample}_R1.fastq.gz",
        r2 = "{fs_prefix}/{dataset}/reads/{query_qc}/{sample}_R2.fastq.gz"
    output:
        bam = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam",
        bai = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam.bai"
    params:
        preset = "default"
    threads: THREADS.get('strobealign', 16)
    wildcard_constraints:
        dataset = "[^/]+",
        reference_path = ".+",
        query_qc = "[^/]+",
        sample = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignment.log"
    benchmark:
        "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/benchmark.txt"
    conda:
        "env.yaml"
    shell:
        """
        echo "=== Strobealign alignment ===" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        echo "Sample: {wildcards.sample}" >> {log}
        echo "Reference path: {wildcards.reference_path}" >> {log}
        echo "Query QC: {wildcards.query_qc}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "Reference: {input.ref}" >> {log}
        echo "Index: {input.index}" >> {log}
        echo "" >> {log}
        
        # Create output directory
        mkdir -p $(dirname {output.bam})
        
        # Align with strobealign, pipe to samtools for sorting
        strobealign \
            -t {threads} \
            {input.ref} \
            {input.r1} \
            {input.r2} \
            2>> {log} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam} \
            - \
            2>> {log}
        
        # Index the BAM
        samtools index -@ {threads} {output.bam} 2>> {log}
        
        echo "" >> {log}
        echo "=== Alignment complete ===" >> {log}
        """


rule strobealign_coverage:
    """
    Calculate per-contig coverage statistics from BAM file.
    Uses samtools coverage for comprehensive metrics.
    """
    input:
        bam = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam",
        bai = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam.bai"
    output:
        coverage = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/coverage.tsv"
    threads: 4
    wildcard_constraints:
        dataset = "[^/]+",
        reference_path = ".+",
        query_qc = "[^/]+",
        sample = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/coverage.log"
    conda:
        "env.yaml"
    shell:
        """
        echo "Calculating coverage" > {log}
        echo "Dataset: {wildcards.dataset}" >> {log}
        echo "Sample: {wildcards.sample}" >> {log}
        
        # Per-contig coverage summary
        samtools coverage {input.bam} > {output.coverage} 2>> {log}
        
        echo "Coverage calculation complete!" >> {log}
        """


rule strobealign_flagstat:
    """
    Generate alignment statistics with samtools flagstat.
    Provides read mapping rates and quality metrics.
    """
    input:
        bam = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam"
    output:
        flagstat = "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/flagstat.txt"
    threads: 1
    wildcard_constraints:
        dataset = "[^/]+",
        reference_path = ".+",
        query_qc = "[^/]+",
        sample = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/flagstat.log"
    conda:
        "env.yaml"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        """
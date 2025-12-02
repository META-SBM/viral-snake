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
    """Find reference file with any common FASTA extension."""
    extensions = ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz']
    for ext in extensions:
        fasta = f"{ref_base}{ext}"
        if os.path.exists(fasta):
            return fasta, ext
    return None, None

def get_reference_fasta(wildcards):
    """
    Construct path to reference FASTA from reference_path wildcard.
    Handles .sanitized suffix automatically.
    """
    ref_path = wildcards.reference_path
    
    # For references/, handle sanitization and extension detection
    if ref_path.startswith("references/"):
        # Check if sanitized version is requested
        if ".sanitized" in ref_path:
            # Extract base path (remove .sanitized)
            ref_base = ref_path.replace(".sanitized", "")
            original, ext = find_reference_file(ref_base)
            if not original:
                raise FileNotFoundError(f"No FASTA file found for {ref_base}")
            return f"{ref_base}.sanitized{ext}"
        else:
            # Return original file
            original, ext = find_reference_file(ref_path)
            if not original:
                raise FileNotFoundError(f"No FASTA file found for {ref_path}")
            return original
    
    # For structured outputs, append standard filenames
    elif ref_path.startswith(("assembly/", "co_assembly/")):
        fasta = f"{ref_path}/contigs.fa"
    elif ref_path.startswith("refseq/"):
        fasta = f"{ref_path}/genomic.fna"
    else:
        fasta = f"{ref_path}/contigs.fa"
    
    # Validation
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Reference FASTA not found: {fasta}")
    
    return fasta


rule sanitize_reference_headers:
    """
    Replace spaces in FASTA headers with ___.
    Works with any FASTA extension (.fa, .fasta, .fna, .gz variants).
    """
    input:
        ref = "references/{ref_path}.{ext}"
    output:
        sanitized = "references/{ref_path}.sanitized.{ext}"
    wildcard_constraints:
        ext = "(fasta|fa|fna)(\.gz)?"
    log:
        "references/{ref_path}.sanitized.{ext}.log"
    conda:
        "env.yaml"
    shell:
        """
        echo "Sanitizing headers in {input.ref}" > {log}
        
        # Handle compressed and uncompressed files
        if [[ {input.ref} == *.gz ]]; then
            zcat {input.ref} | sed '/^>/s/ /___/g' | gzip > {output.sanitized}
        else
            sed '/^>/s/ /___/g' {input.ref} > {output.sanitized}
        fi
        
        echo "Created: {output.sanitized}" >> {log}
        """

def get_strobealign_index(wildcards):
    """Get path to strobealign index file."""
    return f"alignment/__index__/strobealign/{wildcards.reference_path}/r{STROBEALIGN_READ_LENGTH}.sti"

rule strobealign_index:
    """
    Create strobealign index for a reference.
    The index is read-length optimized.
    """
    input:
        ref = get_reference_fasta
    output:
        index = "alignment/__index__/strobealign/{reference_path}/r{read_len}.sti"
    params:
        read_len = lambda wildcards: wildcards.read_len
    threads: 1
    log:
        "alignment/__index__/strobealign/{reference_path}/r{read_len}.log"
    conda:
        "env.yaml"
    shell:
        """
        echo "Creating strobealign index for {wildcards.reference_path}" > {log}
        echo "Read length: {params.read_len}" >> {log}
        echo "Reference: {input.ref}" >> {log}
        
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
        r1 = "reads/{query_qc}/{sample}_R1.fastq.gz",
        r2 = "reads/{query_qc}/{sample}_R2.fastq.gz"
    output:
        bam = "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam",
        bai = "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam.bai"
    params:
        preset = "default"  # Could be expanded for different presets
    threads: THREADS.get('strobealign', 16)
    wildcard_constraints:
        reference_path = ".+",
        query_qc = "[^/]+",
        sample = "[^/]+"
    log:
        "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignment.log"
    benchmark:
        "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/benchmark.txt"
    conda:
        "env.yaml"
    shell:
        """
        echo "Aligning {wildcards.sample} to {wildcards.reference_path}" > {log}
        echo "Query QC: {wildcards.query_qc}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "Reference: {input.ref}" >> {log}
        echo "Index: {input.index}" >> {log}
        
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
        
        echo "Alignment complete!" >> {log}
        """

rule strobealign_coverage:
    """
    Calculate per-contig coverage statistics from BAM file.
    """
    input:
        bam = "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam",
        bai = "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam.bai"
    output:
        coverage = "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/coverage.tsv"
    threads: 4
    log:
        "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/coverage.log"
    conda:
        "env.yaml"
    shell:
        """
        echo "Calculating coverage for {wildcards.sample}" > {log}
        
        # Per-contig coverage summary
        samtools coverage {input.bam} > {output.coverage} 2>> {log}
        
        echo "Coverage calculation complete!" >> {log}
        """

rule strobealign_flagstat:
    """
    Generate alignment statistics with samtools flagstat.
    """
    input:
        bam = "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/alignments.sorted.bam"
    output:
        flagstat = "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/flagstat.txt"
    threads: 1
    log:
        "alignment/strobealign__default/{reference_path}/__reads__/{query_qc}/{sample}/flagstat.log"
    conda:
        "env.yaml"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        """


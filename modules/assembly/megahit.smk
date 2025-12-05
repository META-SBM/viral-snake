rule megahit:
    """
    Assemble reads using MEGAHIT for individual samples.
    Uses temporary directory to avoid cluttering output with intermediate files.
    """
    input:
        r1 = "{fs_prefix}/{dataset}/reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "{fs_prefix}/{dataset}/reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        contigs = "{fs_prefix}/{dataset}/assembly/megahit/{qc_filter}/{sample}/contigs.fa"
    benchmark: 
        "{fs_prefix}/{dataset}/assembly/megahit/{qc_filter}/{sample}/benchmark.txt"
    params:
        basedir = "{fs_prefix}/{dataset}/assembly/megahit/{qc_filter}/{sample}",
        tmpdir = "{fs_prefix}/{dataset}/assembly/megahit/{qc_filter}/{sample}/tmp_megahit",
        intermediate_contigs_dir = "{fs_prefix}/{dataset}/assembly/megahit/{qc_filter}/{sample}/intermediate_contigs"
    threads: THREADS['megahit']
    wildcard_constraints:
        dataset = "[^/]+",
        sample = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/assembly/megahit/{qc_filter}/{sample}/megahit.log"
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        # Create base directory and log
        mkdir -p {params.basedir}
        echo "Assembling {wildcards.dataset}/{wildcards.sample}" > {log}
        
        # Clean up any old temp directory
        rm -rf {params.tmpdir}
        
        # Run MEGAHIT in temp directory
        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.tmpdir} \
            -t {threads} \
            2>> {log}
        
        # Move all MEGAHIT output up one level
        mv {params.tmpdir}/* {params.basedir}/
        
        # Rename final.contigs.fa to contigs.fa
        mv {params.basedir}/final.contigs.fa {output.contigs}
        
        # Clean up empty temp directory
        rm -rf {params.tmpdir}

        # Clean up intermediate_contigs directory
        rm -rf {params.intermediate_contigs_dir}
        """


def get_coassembly_input_files(wildcards):
    """
    Get input read files for co-assembly from collection metadata.
    
    Args:
        wildcards: Must contain fs_prefix, dataset, and collection
        
    Returns:
        dict: Keys 'r1' and 'r2' with lists of read file paths
    """
    import yaml
    
    # Load collection metadata
    collection_file = f"{wildcards.fs_prefix}/{wildcards.dataset}/config/sample_collections/{wildcards.collection}.yaml"
    
    with open(collection_file) as f:
        collection = yaml.safe_load(f)
    
    qc_filter = collection['qc_filter']
    samples = collection['samples']
    
    # Build paths to read files
    r1_files = expand(
        "{fs_prefix}/{dataset}/reads/{qc_filter}/{sample}_R1.fastq.gz",
        fs_prefix=wildcards.fs_prefix,
        dataset=wildcards.dataset,
        qc_filter=qc_filter,
        sample=samples
    )
    
    r2_files = expand(
        "{fs_prefix}/{dataset}/reads/{qc_filter}/{sample}_R2.fastq.gz",
        fs_prefix=wildcards.fs_prefix,
        dataset=wildcards.dataset,
        qc_filter=qc_filter,
        sample=samples
    )
    
    return {'r1': r1_files, 'r2': r2_files}


rule megahit_coassembly:
    """
    Co-assembly using MEGAHIT - combines reads from multiple samples.
    MEGAHIT handles multiple input files natively via comma-separated lists.
    """
    input:
        unpack(get_coassembly_input_files)
    output:
        contigs = "{fs_prefix}/{dataset}/co_assembly/megahit/{collection}/contigs.fa"
    benchmark: 
        "{fs_prefix}/{dataset}/co_assembly/megahit/{collection}/benchmark.txt"
    params:
        basedir = "{fs_prefix}/{dataset}/co_assembly/megahit/{collection}",
        tmpdir = "{fs_prefix}/{dataset}/co_assembly/megahit/{collection}/tmp_megahit",
        r1_list = lambda wildcards, input: ",".join(input.r1),
        r2_list = lambda wildcards, input: ",".join(input.r2)
    threads: 32
    wildcard_constraints:
        dataset = "[^/]+",
        collection = "[^/]+"
    log:
        "{fs_prefix}/{dataset}/co_assembly/megahit/{collection}/megahit.log"
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        # Create base directory and log
        mkdir -p {params.basedir}
        echo "Co-assembling {wildcards.dataset}/{wildcards.collection}" > {log}
        echo "R1 files: {input.r1}" >> {log}
        echo "R2 files: {input.r2}" >> {log}
        
        # Clean up any old temp directory
        rm -rf {params.tmpdir}
        
        # Run MEGAHIT in temp directory
        megahit \
            -1 {params.r1_list} \
            -2 {params.r2_list} \
            -o {params.tmpdir} \
            -t {threads} \
            2>> {log}
        
        # Move all MEGAHIT output up one level
        mv {params.tmpdir}/* {params.basedir}/
        
        # Rename final.contigs.fa to contigs.fa
        mv {params.basedir}/final.contigs.fa {output.contigs}
        
        # Clean up empty temp directory
        rm -rf {params.tmpdir}
        """
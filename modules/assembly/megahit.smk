rule megahit:
    input:
        r1 = "reads/{qc_filter}/{sample}_R1.fastq.gz",
        r2 = "reads/{qc_filter}/{sample}_R2.fastq.gz"
    output:
        contigs = "assembly/megahit/{qc_filter}/{sample}/contigs.fa"
    benchmark: 
        "assembly/megahit/{qc_filter}/{sample}/benchmark.txt"
    params:
        basedir = "assembly/megahit/{qc_filter}/{sample}",
        tmpdir = "assembly/megahit/{qc_filter}/{sample}/tmp_megahit"
    threads: THREADS['megahit']
    log:
        "assembly/megahit/{qc_filter}/{sample}/megahit.log"
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        # Create base directory and log
        mkdir -p {params.basedir}
        echo "Assembling {wildcards.sample}" > {log}
        
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
        """

rule megahit_coassembly:
    """Co-assembly - MEGAHIT handles multiple files natively"""
    input:
        r1 = lambda wildcards: expand(
            "reads/{qc_filter}/{sample}_R1.fastq.gz",
            qc_filter=COLLECTIONS[wildcards.collection]['qc_filter'],
            sample=COLLECTIONS[wildcards.collection]['samples']
        ),
        r2 = lambda wildcards: expand(
            "reads/{qc_filter}/{sample}_R2.fastq.gz",
            qc_filter=COLLECTIONS[wildcards.collection]['qc_filter'],
            sample=COLLECTIONS[wildcards.collection]['samples']
        )
    output:
        contigs = "co_assembly/megahit/{collection}/contigs.fa"
    benchmark: "co_assembly/megahit/{collection}/benchmark.txt"
    params:
        basedir = "co_assembly/megahit/{collection}",
        tmpdir = "co_assembly/megahit/{collection}/tmp_megahit",
        r1_list = lambda wildcards, input: ",".join(input.r1),
        r2_list = lambda wildcards, input: ",".join(input.r2)
    threads: 32
    wildcard_constraints:
        collection = "[^/]+"  # Explicit constraint: no slashes in collection name
    log:
        "co_assembly/megahit/{collection}/megahit.log"
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        # Create base directory and log
        echo "Co-assembling {wildcards.collection}" > {log}
        echo "Samples: {input.r1}" >> {log}
        
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

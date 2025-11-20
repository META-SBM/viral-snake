rule mmseqs_deduplicate:
    input:
        contigs = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/contigs.fa"
    output:
        rep_seqs = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/deduplicated_id{identity}_cov{coverage}/contigs_deduplicated.fa",
        cluster_tsv = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/deduplicated_id{identity}_cov{coverage}/cluster_info.tsv",
        all_seqs = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/deduplicated_id{identity}_cov{coverage}/all_seqs.fa"
    params:
        prefix = "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/deduplicated_id{identity}_cov{coverage}/result",
        min_seq_id = lambda wildcards: float(wildcards.identity) / 100,
        coverage = lambda wildcards: float(wildcards.coverage) / 100,
        cov_mode = 2,
        cluster_mode = 0,
        sensitivity = 4.0
    threads: 16
    log:
        "assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/deduplicated_id{identity}_cov{coverage}/log"
    conda:
        "mmseqs2"
    shell:
        """
        TMPDIR=$(mktemp -d -p . mmseqs_tmp.XXXXXX)
        
        mmseqs easy-cluster {input.contigs} {params.prefix} $TMPDIR \
            --min-seq-id {params.min_seq_id} \
            -c {params.coverage} \
            --cov-mode {params.cov_mode} \
            --cluster-mode {params.cluster_mode} \
            -s {params.sensitivity} \
            --threads {threads} \
            -v 3 \
            2>&1 | tee {log}
        
        mv {params.prefix}_rep_seq.fasta {output.rep_seqs}
        mv {params.prefix}_cluster.tsv {output.cluster_tsv}
        mv {params.prefix}_all_seqs.fasta {output.all_seqs}
        
        rm -rf $TMPDIR
        """
# Snakefile for tick metavirome assembly pipeline

import yaml
import glob
import os

# Including all submodules
include: "modules/mmseqs2/cluster.smk"
include: "modules/assembly/anvio.smk"
include: "modules/diamond/diamond.smk"
include: "modules/minimap2.smk"

include: "modules/qc/cutadapt.smk"
include: "modules/qc/fastqc.smk"
include: "modules/qc/read_stats.smk"

include: "modules/kraken2/kraken2.smk"
include: "modules/kraken2/kraken2_contigs.smk"
include: "modules/blast/blast.smk"

include: "modules/assembly/megahit.smk"
include: "modules/assembly/metaspades.smk"
include: "modules/assembly/qc.smk"


# Get workdir from config or use default
WORKDIR = config.get("workdir")

workdir: WORKDIR

# SAMPLES = ["KLP_168"]
# print(os.path.join(WORKDIR, 'reads/raw', '*_1.fq'))
# SAMPLES = glob.glob(os.path.join(WORKDIR, 'reads/raw', '*_1.fq'))
# SAMPLES = [s.split('/')[-1].replace('_1.fq', '') for s in SAMPLES]

# print(os.path.join(WORKDIR, 'reads/raw', '*_1.fq'))
SAMPLES = glob.glob(os.path.join(WORKDIR, 'reads/raw', '*_R1.fastq.gz'))
SAMPLES = [s.split('/')[-1].replace('_R1.fastq.gz', '') for s in SAMPLES]

# print(SAMPLES)
QC_FILTER = "raw__cutadapt_no_mgi_min_len_90"
ASSEMBLERS = ["megahit", "metaspades", "megahit_rna"]
ASSEMBLERS = ["megahit"]


# Discover collections
COLLECTION_DIR = Path(WORKDIR) / "sample_collections"
COLLECTIONS = {}
if COLLECTION_DIR.exists():
    for yaml_file in COLLECTION_DIR.glob("*.yaml"):
        with open(yaml_file) as f:
            coll = yaml.safe_load(f)
            COLLECTIONS[coll['name']] = coll

COOLS = ['gus_khrustalny_district__novoopokino_village_dermacentor_reticulatus', 'gvardeysky_district__konstantinovka_settlement_dermacentor_reticulatus',	'chuguyevsky_district__zhuravlevka_river_ixodes_persulcatus', 'chuguyevsky_district__zhuravlevka_river_haemophysalis_japonica', 'blagoveshchensky_district__jsc_polief_dermacentor']

refs = glob.glob(os.path.join(WORKDIR, 'refseq_reference/*', 'ncbi_dataset/data/genomic.fna'))
refs = [r.split('/')[-4] for r in refs]
# print(refs)


rule all:
    input:
        # FastQC reports
        expand("qc/read_stats/{qc_filter}/{sample}_read_counts.tsv",
               qc_filter=QC_FILTER, sample=SAMPLES),
        # Co-assemblies
        # expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/contigs.fa",
        #        assembler=ASSEMBLERS, collection=list(COLLECTIONS.keys()), min_len=3000),
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/diamond_faster/NR/hits_with_taxonomy.tsv",
               assembler=ASSEMBLERS, collection=COLLECTIONS, min_len=800),
        # expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/contig_stats.tsv",
        #        assembler=ASSEMBLERS, collection=list(COLLECTIONS.keys()), min_len=800),
        expand("co_assembly/{assembler}/{collection}/contigs_formatted_minlen_{min_len}/minimap2/{reference_org}/alignments.sorted.bam",
               assembler=ASSEMBLERS, collection=COLLECTIONS, min_len=800, reference_org=refs)

        # Assemblies
        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs.fa",
        #        assembler=ASSEMBLERS, qc_filter=QC_FILTER, sample=SAMPLES),
        # Individual QUAST reports
        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/contig_stats.tsv",
        #        assembler=ASSEMBLERS, qc_filter=QC_FILTER, sample=SAMPLES, min_len=300),
        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/kraken2/confidence_{confidence}/kraken2_output_with_taxonomy.tsv",
        #        assembler=ASSEMBLERS, qc_filter=QC_FILTER, sample=SAMPLES, min_len=300,confidence=0.01),
        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/diamond/{database}/hits_with_taxonomy.tsv",
        #        assembler=ASSEMBLERS, qc_filter=QC_FILTER, sample=SAMPLES[0], min_len=1000, database="NR"),

        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/blast/{database}/blastn.tsv",
        #        assembler=ASSEMBLERS, qc_filter=QC_FILTER, sample=SAMPLES, min_len=300, database="core_nt"),
        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/blast/{database}/blastn_with_taxonomy.tsv",
        #        assembler=ASSEMBLERS, qc_filter=QC_FILTER, sample=SAMPLES[0], min_len=300, database="core_nt"),



        # expand("kraken2/{qc_filter}/{sample}.bracken.report",
        #        qc_filter=QC_FILTER, sample=SAMPLES),
        # "feature_tables/bracken-species-all-0.4-min-len-90/taxonomy_table.tsv",
        # Comparison QUAST report
        # expand("qc/quast/comparison/{qc_filter}/{sample}/report.html",
        #        qc_filter=QC_FILTER, sample=SAMPLES),
    #     expand("annotation/virsorter2/megahit/{qc_filter}/{sample}/final-viral-combined.fa",
	#    qc_filter=QC_FILTER, sample=SAMPLES)

        # expand("assembly/{assembler}/{qc_filter}/{sample}/contigs_formatted_minlen_{min_len}/deduplicated_id95_cov85/contigs_deduplicated.fa",
        #     assembler=ASSEMBLERS, qc_filter=QC_FILTER, sample=SAMPLES[0], min_len=300)

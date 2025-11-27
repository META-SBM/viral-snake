import glob
import os

# Base paths
BASE_DIR = "/mnt/mgx/DATASETS/INTERNAL/VIROME/VIROME4/F350035952"
L01_DIR = os.path.join(BASE_DIR, "L01")
L02_DIR = os.path.join(BASE_DIR, "L02")
MERGED_DIR = os.path.join(BASE_DIR, "MERGED_LANES")

# Auto-detect sample IDs from L01 files
L01_files = glob.glob(os.path.join(L01_DIR, "F350035952_L01_*_1.fq.gz"))
SAMPLES = [os.path.basename(f).replace("F350035952_L01_", "").replace("_1.fq.gz", "") for f in L01_files]

rule all:
    input:
        expand(os.path.join(MERGED_DIR, "F350035952_{sample}_{read}.fq.gz"), 
               sample=SAMPLES, read=[1, 2])

rule merge_lanes:
    input:
        l01 = os.path.join(L01_DIR, "F350035952_L01_{sample}_{read}.fq.gz"),
        l02 = os.path.join(L02_DIR, "F350035952_L02_{sample}_{read}.fq.gz")
    output:
        os.path.join(MERGED_DIR, "F350035952_{sample}_{read}.fq.gz")
    shell:
        "cat {input.l01} {input.l02} > {output}"
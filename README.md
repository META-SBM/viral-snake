# viral-snake
Snakemake pipelines for metavirome data processing, visualisation and analysis


# Установка

## Установка пакета
`pip install -e ./`

## Потребуется утсановить anvio-8

```
conda create -y --name anvio-8 python=3.10

conda activate anvio-8

conda install -y -c conda-forge -c bioconda python=3.10 \
        sqlite=3.46 prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript \
        nodejs=20.12.2

conda install -y -c bioconda fastani

curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz \
        --output anvio-8.tar.gz

pip install anvio-8.tar.gz
```

# Запуск

Нужно указать флаг `--use-conda`
И проставить `--conda-prefix` - тогда окружения будут жить в этой папочке и не пересоздаваться при смене рабочей директории

`snakemake -s Snakefile.smk --use-conda --cores 320 --config workdir=/mnt/mgx/DATASETS/INTERNAL/VIROME/TEST --rerun-incomplete --latency-wait 100 --conda-prefix /mnt/mgx/RUNS/fedorov_de/SNAKE-CONDA --rerun-triggers mtime`


# Менеджмент

`viral-snake info samples /mnt/mgx/DATASETS/INTERNAL/VIROME/VIROME4`
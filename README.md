# viral-snake
Snakemake pipelines for metavirome data processing, visualisation and analysis




# Запуск

Нужно указать флаг `--use-conda`
И проставить `--conda-prefix` - тогда окружения будут жить в этой папочке и не пересоздаваться

`snakemake -s Snakefile.smk --use-conda --cores 320 --config workdir=/mnt/mgx/DATASETS/INTERNAL/VIROME/RUN3 --rerun-incomplete --latency-wait 100 --conda-prefix /mnt/mgx/DATASETS/INTERNAL/VIROME/RUN3`
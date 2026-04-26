#!/bin/bash
# Activate before running: conda activate fastqc_env

# Download SRA run using SRAtools
fasterq-dump -e 1 -p SRR36691286

# Compress into a gzip
gzip SRR36691286_1.fastq
gzip SRR36691286_2.fastq

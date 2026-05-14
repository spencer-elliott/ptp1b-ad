#!/bin/bash
# Activate before running: conda activate fastqc_env

# Make a new directory for data, and another inside for fastq files
mkdir -p data/fastq

# Download SRA run using SRAtools
fasterq-dump -e 1 -p SRR36691286 -O data/fastq

# Compress into a gzip
gzip data/fastq/SRR36691286_1.fastq
gzip data/fastq/SRR36691286_2.fastq

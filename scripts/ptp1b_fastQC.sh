#!/bin/bash
# Activate before running: conda activate fastqc_env

# Create an output directory for fastqc
mkdir -p data/fastqc

# Run fastqc to see the quality of raw data/reads
fastqc data/fastq/SRR36691286_1.fastq.gz data/fastq/SRR36691286_2.fastq.gz -o data/fastqc

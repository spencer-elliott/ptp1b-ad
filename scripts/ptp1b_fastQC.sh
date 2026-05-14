#!/bin/bash
# Activate before running: conda activate fastqc_env

# Run fastqc to see the quality of raw data/reads
fastqc SRR36691286_1.fastq.gz SRR36691286_2.fastq.gz

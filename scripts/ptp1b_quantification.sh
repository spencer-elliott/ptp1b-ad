#!/bin/bash
# Activate before running: conda activate kallisto_env
# Files needed: Reference transcriptome, 10xv3 whitelist, transcripts.txt and transcripts_to_genes.txt 
#   See corresponding Jupyter notebook for instructions

# Create directories for reference files and kallisto output
mkdir -p data/references
mkdir -p results/kallisto_output

# Create an index for kallisto from the mouse reference transcriptome:
kallisto index -i data/references/mm10.idx \
  data/references/gencode.vM23.transcripts.fa.gz

# Pseudoalignment of reads (for 10x data)
kallisto bus \
  -i data/references/mm10.idx \
  -o results/kallisto_output \
  -x 10xv3 \
  -t 4 \
  data/fastq/SRR36691286_1.fastq.gz data/fastq/SRR36691286_2.fastq.gz

# Ordering/sorting the barcodes
bustools sort -t 4 \
  -o results/kallisto_output/sorted.bus \
  results/kallisto_output/output.bus

# Fixing sequencing errors (needs whitelist)
bustools correct \
  -w data/references/3M-february-2018.txt \
  -o results/kallisto_output/corrected.bus \
  results/kallisto_output/sorted.bus

# Re-sorting after correction
bustools sort \
  -t 4 \
  -o results/kallisto_output/corrected_sorted.bus \
  results/kallisto_output/corrected.bus

# Creating a gene expression matrix - gives us gene expression counts per cell
bustools count \
  -o results/kallisto_output/counts \
  -g results/kallisto_output/transcripts_to_genes.txt \
  -e results/kallisto_output/matrix.ec \
  -t results/kallisto_output/transcripts.txt \
  --genecounts \
  results/kallisto_output/corrected_sorted.bus

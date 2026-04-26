#!/bin/bash
# Activate before running: conda activate kallisto_env

# Create an index for kallisto from the mouse reference transcriptome:
kallisto index -i mm10.idx gencode.vM23.transcripts.fa.gz

# Pseudoalignment of reads (for 10x data)
kallisto bus -i mm10.idx -o kallisto_output -x 10xv3 -t 4 SRR36691286_1.fastq.gz SRR36691286_2.fastq.gz

# Ordering/sorting the barcodes
bustools sort -t 4 -o kallisto_output/sorted.bus kallisto_output/output.bus

# Fixing sequencing errors (needs whitelist)
bustools correct -w 3M-february-2018.txt -o /dev/stdout kallisto_output/sorted.bus | \
  bustools sort -t 4 -o kallisto_output/corrected_sorted.bus -

# Creating a gene expression matrix - gives us gene expression counts per cell
bustools count \
  -o kallisto_output/counts \
  -g kallisto_output/transcripts_to_genes.txt \
  -e kallisto_output/matrix.ec \
  -t kallisto_output/transcripts.txt \
  kallisto_output/corrected.bus

bustools count \
  -o kallisto_output/counts \
  -g kallisto_output/transcripts_to_genes.txt \
  -e kallisto_output/matrix.ec \
  -t kallisto_output/transcripts.txt \
  --genecounts \
  kallisto_output/corrected_sorted.bus

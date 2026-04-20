#!/bin/bash

# Create an index for kallisto from the mouse reference transcriptome:
kallisto index -i mm10.idx gencode.vM38.transcripts.fa.gz

# Pseudoalignment of reads (for 10x data)
kallisto bus -i mm10.idx -o kallisto_output -x 10xv3 -t 4 SRR36691286_1.fastq.gz SRR36691286_2.fastq.gz

# ON THIS STEP 4/20/26!
# Ordering/sorting the barcodes
bustools sort -t 4 -o kallisto_output/sorted.bus kallisto_output/output.bus

# Fixing sequencing errors
bustools correct -w whitelist.txt -o kallisto_output/corrected.bus kallisto_output/sorted.bus

# Creating a gene expression matrix
bustools count -o kallisto_output/counts \
  -g transcripts_to_genes.txt \
  -e matrix.ec \
  -t transcripts.txt \
  kallisto_output/corrected.bus

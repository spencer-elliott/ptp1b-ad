# RNA-Seq Pipeline

This is my own practice playing around with data from single cell RNA-seq. This data was in a paper that Yuxin from my lab released last year:

*"PTP1B inhibition promotes microglial phagocytosis  in Alzheimer’s disease models by enhancing SYK signaling"* [link]

My goal is to understand and recreate as much analysis as I can! And maybe check out some other cell populations, statistics as well. It's a good chance to try out some awesome tools that are used on this type of data.

## Background

*How does PTP1B inhibition affect microglial gene expression in APP/PS1 mice?*

I wrote a summary of this research question that you can find in BACKGROUND.md!

This scRNA-seq data is from the brains of two groups of APP/PS1 mice, a widely used transgenic model for Alzheimer's disease:
- APP/PS1-PTP1B+/+ (Alzheimer's model with WT PTP1B)
- APP/PS1-PTP1B-/- (Alzheimer's model with PTP1B knockout)

Four mice were chosen from each experimental group. Each sample for RNA sequencing was pooled together from two mice, for a total of four samples.

## Pipeline

1. 'ptp1b_download.sh' - Download raw FASTQ data using SRA Toolkit
2. 'ptp1b_fastQC.sh' - Quality control using FastQC
3. 'ptp1b_quantification.sh' - Pseudoalignment using kallisto/bustools

Alongside each of these scripts I'm going to create a companion Jupyter notebook that documents each choice I made, and what each of these commands is doing.

## Stuff I learned

- Conda
- Kallisto/bustools
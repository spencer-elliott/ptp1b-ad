# RNA-Seq Pipeline

This is my own practice playing around with data from single cell RNA-seq. This data is from a paper that Yuxin from my lab released last year:

*"PTP1B inhibition promotes microglial phagocytosis  in Alzheimer’s disease models by enhancing SYK signaling"* [link]

My goal is to understand/recreate as much analysis as I can! And maybe check out some other cell populations, statistics as well. And try out some of the awesome tools that are used for this type of data.

## Background

*How does PTP1B inhibition affect microglial gene expression in APP/PS1 mice?*

I wrote a summary of this reserch question that you can find in BACKGROUND.md. scRNA-seq was performed on the brains of two populations of APP/PS1 mice, a widely used transgenic model of Alzheimer's disease:
- APP/PS1-PTP1B+/+ (Alzheimer's model with WT PTP1B)
- APP/PS1-PTP1B-/- (Alzheimer's model with PTP1B knockout)

Each of these groups contains 4 mice, and each sample for RNA sequencing was pooled from 2 mice.

The main cell type we are investigating is the microglia, and the effects of PTP1B inhibition on microglial gene expression.

## Methods

You can see a full description of my thought process in METHODS.md.

## Stuff I learned

- Conda
- Kallisto/bustools
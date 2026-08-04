# RNA-Seq Pipeline

This is my own practice playing around with data from single cell RNA-seq. The data was collected for a paper that Yuxin from my lab released last year:

*"PTP1B inhibition promotes microglial phagocytosis in Alzheimer’s disease models by enhancing SYK signaling"* [link]

scRNA-seq was performed to investigate gene expression across various cell types. My goal is to understand and recreate as much of their analysis pipeline as I can. It's a good chance to try the neat tools that are used on this kind of data!

## Background

*How does PTP1B inhibition affect microglial gene expression in APP/PS1 mice?*

One of the main symptoms of Alzheimer's disease is an accumulation of Amyloid β proteins into plaques that are toxic to the brain. Microglia are a type of immune cell found in the brain that can remove these plaques via phagocytosis, making them a key target for treatment research.

Protein-tyrosine phosphatase 1B (PTP1B) is an enzyme with several functions in metabolic control, as well as proven roles in immune cell signaling. Studies linking metabolic regulation to microglial function, as well as the loss of PTP1B to increased B cell activation, suggested the potential for PTP1B to have a role in microglial function, specifically the SYK pathway in microglia that affects both metabolic fitness and phagocytosis.

In this paper PTP1B was found to be highly expressed in microglia and have a role directly acting on the SYK pathway. Inhibiting PTP1B was found to promote phagocytosis and the immune response targeting Amyloid β plaques, reducing amyloid levels and improving cognitive function in mouse models of Alzheimer's disease.

This scRNA-seq data is from the brains of two groups of APP/PS1 mice, a widely used transgenic model for Alzheimer's disease:
- APP/PS1-PTP1B+/+ (Alzheimer's model with WT PTP1B)
- APP/PS1-PTP1B-/- (Alzheimer's model with PTP1B knockout)

Four mice were chosen from each experimental group. Each sample for RNA sequencing was pooled together from two mice, for a total of four samples.

## Pipeline/Steps

For each of these steps I have a script as well as a corresponding Jupyter notebook (for shell scripts) or rendered Quarto file (for R analysis) that documents each choice I made.

1. Downloading raw FASTQ data using SRA Toolkit
- Script: ['ptp1b_download.sh'](scripts/ptp1b_download.sh)
- Notebook: ['ptp1b_download.ipynb'](notebooks/ptp1b_download.ipynb)

2. Quality control using FastQC
- Script: ['ptp1b_fastQC.sh'](scripts/ptp1b_fastQC.sh)
- Notebook: ['ptp1b_fastQC.ipynb'](notebooks/ptp1b_fastQC.ipynb)

3. Pseudoalignment and quantification for one sample using Kallisto/Bustools
- Script: ['ptp1b_quantification.sh'](scripts/ptp1b_quantification.sh)
- Notebook: ['ptp1b_quantification.ipynb'](notebooks/ptp1b_quantification.ipynb)

4. Analysis of a feature-barcode matrix in R using Seurat
- Script: ['ptp1b_seurat.R'](scripts/ptp1b_seurat.R)
- Analysis: ['ptp1b_seurat.html'](https://spencerelliott.github.io/ptp1b-ad/notebooks/ptp1b_seurat.html)

5. Alignment and quantification for all four samples using CellRanger
- Scripts: ['slurm' folder](scripts/slurm/)
- Notebook: ['ptp1b_slurm.ipynb'](notebooks/ptp1b_slurm.ipynb)

6. 'ptp1b_seurat_final.R' - Analysis for all four samples in R using Seurat
- Script: ['ptp1b_seurat_final.R'](scripts/ptp1b_seurat_final.R)
- Analysis: ['ptp1b_seurat_final.html'](https://spencerelliott.github.io/ptp1b-ad/notebooks/ptp1b_seurat_final.html)

## Stuff I used

- Conda
- Kallisto
- bustools
- CellRanger
- CellBender
- HPC and Slurm task files
- Seurat
- SCTransform
- Harmony
- SingleR
- Manual annotation using marker genes

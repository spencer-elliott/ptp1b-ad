#!/bin/bash
#SBATCH --job-name=ptp1b_cellbender
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -e /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_cellbender.err
#SBATCH -o /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_cellbender.out

module load EBModules
module load Anaconda3

source /grid/it/data/elzar/easybuild/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cellbender

cd ~/ptp1b_ad
mkdir -p results/cellbender_output/KO_SRR36691286
mkdir -p results/cellbender_output/KO_SRR36691287
mkdir -p results/cellbender_output/AD_SRR36691288
mkdir -p results/cellbender_output/AD_SRR36691289

cellbender remove-background \
    --input ~/ptp1b_ad/results/cellranger_output/KO_SRR36691286/outs/raw_feature_bc_matrix.h5 \
    --output ~/ptp1b_ad/results/cellbender_output/KO_SRR36691286_cellbender.h5 \
    --epochs 150 \
    --cpu-threads 8 \
    --checkpoint-mins 9999

cellbender remove-background \
    --input ~/ptp1b_ad/results/cellranger_output/KO_SRR36691287/outs/raw_feature_bc_matrix.h5 \
    --output ~/ptp1b_ad/results/cellbender_output/KO_SRR36691287_cellbender.h5 \
    --epochs 150 \
    --cpu-threads 8 \
    --checkpoint-mins 9999

cellbender remove-background \
    --input ~/ptp1b_ad/results/cellranger_output/AD_SRR36691288/outs/raw_feature_bc_matrix.h5 \
    --output ~/ptp1b_ad/results/cellbender_output/AD_SRR36691288_cellbender.h5 \
    --epochs 150 \
    --cpu-threads 8 \
    --checkpoint-mins 9999

cellbender remove-background \
    --input ~/ptp1b_ad/results/cellranger_output/AD_SRR36691289/outs/raw_feature_bc_matrix.h5 \
    --output ~/ptp1b_ad/results/cellbender_output/AD_SRR36691289_cellbender.h5 \
    --epochs 150 \
    --cpu-threads 8 \
    --checkpoint-mins 9999

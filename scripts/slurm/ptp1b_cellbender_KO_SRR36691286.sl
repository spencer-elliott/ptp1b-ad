#!/bin/bash
#SBATCH --job-name=ptp1b_cellbender_KO_SRR36691286
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH -e /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_cellbender_KO_SRR36691286.err
#SBATCH -o /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_cellbender_KO_SRR36691286.out

module load EBModules
module load Anaconda3

source /grid/it/data/elzar/easybuild/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cellbender

cd ~/ptp1b_ad
mkdir -p results/cellbender_output/KO_SRR36691286

cellbender remove-background \
    --input ~/ptp1b_ad/results/cellranger_output/KO_SRR36691286/outs/raw_feature_bc_matrix.h5 \
    --output ~/ptp1b_ad/results/cellbender_output/KO_SRR36691286_cellbender.h5 \
    --epochs 150 \
    --cpu-threads 8 \
    --checkpoint-mins 9999

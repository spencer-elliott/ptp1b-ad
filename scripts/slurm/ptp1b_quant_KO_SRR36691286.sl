#!/bin/bash
#SBATCH --job-name=ptp1b_quant_KO_SRR36691286
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -e /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_quant_KO_SRR36691286.err
#SBATCH -o /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_quant_KO_SRR36691286.out

module load EB5
module load CellRanger/9.0.1

cd ~/ptp1b_ad
mkdir -p results/cellranger_output/KO_SRR36691286

cellranger count \
  --id=KO_SRR36691286 \
  --transcriptome=data/references/refdata-gex-mm10-2020-A \
  --fastqs=data/fastq \
  --sample=KO_SRR36691286 \
  --localcores=8 \
  --localmem=64 \
  --output-dir=results/cellranger_output/KO_SRR36691286

#!/bin/bash
#SBATCH --job-name=ptp1b_quantification
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -e /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_quantification.err
#SBATCH -o /grid/tonks/home/elliott/ptp1b_ad/logs/ptp1b_quantification.out

module load EB5
module load CellRanger/9.0.1

cd ~/ptp1b_ad
mkdir -p results/cellranger_output/KO_SRR36691286
mkdir -p results/cellranger_output/KO_SRR36691287
mkdir -p results/cellranger_output/AD_SRR36691288
mkdir -p results/cellranger_output/AD_SRR36691289

cellranger count \
  --id=KO_SRR36691286 \
  --transcriptome=data/references/refdata-gex-mm10-2020-A \
  --fastqs=data/fastq \
  --sample=KO_SRR36691286 \
  --localcores=8 \
  --localmem=64 \
  --output-dir=results/cellranger_output/KO_SRR36691286 \
  --create-bam false

cellranger count \
  --id=KO_SRR36691287 \
  --transcriptome=data/references/refdata-gex-mm10-2020-A \
  --fastqs=data/fastq \
  --sample=KO_SRR36691287 \
  --localcores=8 \
  --localmem=64 \
  --output-dir=results/cellranger_output/KO_SRR36691287 \
  --create-bam false

cellranger count \
  --id=AD_SRR36691288 \
  --transcriptome=data/references/refdata-gex-mm10-2020-A \
  --fastqs=data/fastq \
  --sample=AD_SRR36691288 \
  --localcores=8 \
  --localmem=64 \
  --output-dir=results/cellranger_output/AD_SRR36691288 \
  --create-bam false

cellranger count \
  --id=AD_SRR36691289 \
  --transcriptome=data/references/refdata-gex-mm10-2020-A \
  --fastqs=data/fastq \
  --sample=AD_SRR36691289 \
  --localcores=8 \
  --localmem=64 \
  --output-dir=results/cellranger_output/AD_SRR36691289 \
  --create-bam false

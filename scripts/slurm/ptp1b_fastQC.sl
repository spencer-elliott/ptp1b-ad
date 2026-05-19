#!/bin/bash
#SBATCH --job-name=ptp1b_fastQC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=10:00:00
#SBATCH -e ~/ptp1b_ad/logs/ptp1b_fastQC.err
#SBATCH -o ~/ptp1b_ad/logs/ptp1b_fastQC.out

module load EB5
module load FastQC

cd ~/ptp1b_ad
mkdir -p data/fastqc

fastqc data/fastq/KO_SRR36691286_1.fastq.gz data/fastq/KO_SRR36691286_2.fastq.gz -o data/fastqc
fastqc data/fastq/KO_SRR36691287_1.fastq.gz data/fastq/KO_SRR36691287_2.fastq.gz -o data/fastqc
fastqc data/fastq/AD_SRR36691288_1.fastq.gz data/fastq/AD_SRR36691288_2.fastq.gz -o data/fastqc
fastqc data/fastq/AD_SRR36691289_1.fastq.gz data/fastq/AD_SRR36691289_2.fastq.gz -o data/fastqc

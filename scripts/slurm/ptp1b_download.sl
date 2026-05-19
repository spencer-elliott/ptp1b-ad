#!/bin/bash
#SBATCH --job-name=ptp1b_download
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=10:00:00
#SBATCH -e ~/ptp1b_ad/logs/ptp1b_download.err
#SBATCH -o ~/ptp1b_ad/logs/ptp1b_download.out

module load EB5
module load SRA-Toolkit
prefetch -O ~/ptp1b_ad/data/sra SRR36691286 SRR36691287 SRR36691288 SRR36691289

cd ~/ptp1b_ad
mkdir -p data/fastq

fasterq-dump -e 4 -p -O data/fastq data/sra/SRR36691286/SRR36691286.sra
mv data/fastq/SRR36691286_1.fastq data/fastq/KO_SRR36691286_1.fastq
mv data/fastq/SRR36691286_2.fastq data/fastq/KO_SRR36691286_2.fastq
gzip data/fastq/KO_SRR36691286_1.fastq
gzip data/fastq/KO_SRR36691286_2.fastq

fasterq-dump -e 4 -p -O data/fastq data/sra/SRR36691287/SRR36691287.sra
mv data/fastq/SRR36691287_1.fastq data/fastq/KO_SRR36691287_1.fastq
mv data/fastq/SRR36691287_2.fastq data/fastq/KO_SRR36691287_2.fastq
gzip data/fastq/KO_SRR36691287_1.fastq
gzip data/fastq/KO_SRR36691287_2.fastq

fasterq-dump -e 4 -p -O data/fastq data/sra/SRR36691288/SRR36691288.sra
mv data/fastq/SRR36691288_1.fastq data/fastq/AD_SRR36691288_1.fastq
mv data/fastq/SRR36691288_2.fastq data/fastq/AD_SRR36691288_2.fastq
gzip data/fastq/AD_SRR36691288_1.fastq
gzip data/fastq/AD_SRR36691288_2.fastq

fasterq-dump -e 4 -p -O data/fastq data/sra/SRR36691289/SRR36691289.sra
mv data/fastq/SRR36691289_1.fastq data/fastq/AD_SRR36691289_1.fastq
mv data/fastq/SRR36691289_2.fastq data/fastq/AD_SRR36691289_2.fastq
gzip data/fastq/AD_SRR36691289_1.fastq
gzip data/fastq/AD_SRR36691289_2.fastq

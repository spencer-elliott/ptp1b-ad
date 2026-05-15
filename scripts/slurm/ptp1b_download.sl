#!/bin/tcsh
#SBATCH --job-name=ptp1b_download
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=10:00:00
#SBATCH -e ptp1b_download.err
#SBATCH -o ptp1b_download.out

cd ptp1b_ad

fasterq-dump -e 1 -p SRR36691286
mv SRR36691286_1.fastq KO_SRR36691286_1.fastq
mv SRR36691286_2.fastq KO_SRR36691286_2.fastq
gzip data/fastq/KO_SRR36691286_1.fastq
gzip data/fastq/KO_SRR36691286_2.fastq

fasterq-dump -e 1 -p SRR36691287
mv SRR36691287_1.fastq KO_SRR36691287_1.fastq
mv SRR36691287_2.fastq KO_SRR36691287_2.fastq
gzip data/fastq/KO_SRR36691287_1.fastq
gzip data/fastq/KO_SRR36691287_2.fastq

fasterq-dump -e 1 -p SRR36691288
mv SRR36691288_1.fastq AD_SRR36691288_1.fastq
mv SRR36691288_2.fastq AD_SRR36691288_2.fastq
gzip data/fastq/AD_SRR36691288_1.fastq
gzip data/fastq/AD_SRR36691288_2.fastq

fasterq-dump -e 1 -p SRR36691289
mv SRR36691289_1.fastq AD_SRR36691289_1.fastq
mv SRR36691289_2.fastq AD_SRR36691289_2.fastq
gzip data/fastq/AD_SRR36691289_1.fastq
gzip data/fastq/AD_SRR36691289_2.fastq

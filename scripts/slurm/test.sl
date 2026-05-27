#!/bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=00:05:00
#SBATCH -e /grid/tonks/home/elliott/ptp1b_ad/logs/test.err
#SBATCH -o /grid/tonks/home/elliott/ptp1b_ad/logs/test.out

echo "job started"
module load EB5
module load CellRanger/9.0.1
echo "modules loaded"
cellranger --version
echo "done"
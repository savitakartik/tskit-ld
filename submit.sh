#! /bin/bash

#SBATCH --job-name=ld_calc_comparisons
#SBATCH --partition=compute
#SBATCH --mem=60G
#SBATCH --array=6-10
#SBATCH --output=logs/slurm/%a.out
#SBATCH --error=logs/slurm/%a.err
#SBATCH --time=300:00:00
#SBACTH --mail-type=END
#SBATCH --mail-user=osvk@novonordisk.com

. /etc/profile.d/modules.sh
module load plink/1.90
conda activate ~/.conda/envs/msprime-tskit-env/

date

i=${SLURM_ARRAY_TASK_ID}
i=$((i-1))
sample_sizes=(10 100 1000 10000 25000 50000 75000 100000 250000 500000)
#sample_sizes=(10 100 1000 10000 25000 50000)
sample_size=${sample_sizes[i]}
echo -ne "Sample size: ${sample_size}\n"

/usr/bin/time -v python tskit_plink_comparisons_batch.py -n ${sample_size}

date

# iMMC

#!/bin/bash
#SBATCH --job-name=sequeez
#SBATCH --output=./Log/sequeez.%A_%a.out
#SBATCH --error=./Log/sequeez.%A_%a.err
#SBATCH --time=05-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=249GB
#SBATCH --array=1-73

module load squeezemeta/1.6.0

cd /ibex/project/c2207/iMMC/Metagenome/squeezeMeta

sample=`cat sample.list | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

cd ${sample}

sqm_longreads.pl -p ${sample} -s test.samples -f data/ -t 32


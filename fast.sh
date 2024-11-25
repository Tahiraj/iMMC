#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --output=./LogF/fastp.%A_%a.out
#SBATCH --error=./LogF/fastp.%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=32
#SBATCH --array=1-73

sample=`ls data/fastq|grep .fastq.gz| head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

prefix=`basename $sample .fastq.gz`

current=$(pwd)
echo $current

mkdir -p "PreProcess/$prefix"

samplePath="data/fastq"

PProcDir="PreProcess/${prefix}";

r1="${samplePath}/${prefix}.fastq.gz"

PorechopR1="${PProcDir}/${prefix}_porechop.fastq.gz";

#Activate conda environment in node terminal as well
conda init bash 
conda activate /home/jamilt/miniconda3/envs/Porechop
porechop -i $r1 -o ${PorechopR1} --threads 32
conda deactivate

module load fastp
filtR1="${PProcDir}/${prefix}.filt.fastq.gz";
echo $r1;   echo $filtR1; 

fastp -i ${PorechopR1} -o ${filtR1} -q 12 --thread 32 --dedup -j ${PProcDir}/${prefix}_fastp.json -h ${PProcDir}/${prefix}_fastp.html;

module load fastqc/0.12.0
mkdir -p $current/PreProcess_q12/qscores
qscore="$current/PreProcess_q12/qscores"
time -p fastqc --threads 32 ${filtR1} -o ${qscore}

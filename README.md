# iMMC
Exploring microbial diversity using portable technology for in situ analysis.

For metagenomic analysis we have 72 samples of Coral and Mangroves and one positive control. 
The sample file specifies the samples, the names of their corresponding raw read files and the sequencing pair represented in those files, separated by tabulators.
It has the format: 
```
<Sample>	 <filename> 	<pair1|pair2>	
```
Prepare one sample file  
```
C1	C1.filt.fastq.gz	pair1
```
and then replicate for each sample from a list by replacing the sample name  with new sample name.
For single long read use only `pair1`  
The code below create a new sample file and also replaces the `C1` with `C2`.

```
#!/bin/bash
#SBATCH --job-name=sequeez
#SBATCH --output=./Log/sequeez.%A_%a.out
#SBATCH --error=./Log/sequeez.%A_%a.err
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --array=1-73

cd /ibex/project/c2207/iMMC/Metagenome

sample=`cat sample.list | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

prefix=`basename $sample .fastq.gz`

mkdir -p  squeezeMeta/${prefix}/data

cp PreProcess_q12/${prefix}/${prefix}.filt.fastq.gz squeezeMeta/${prefix}/data

cp SqueezeMeta/test.samples squeezeMeta/${prefix}/

cd squeezeMeta/${prefix}

oldstring=“C1"

newstring=“$prefix"

sed -i "s/$oldstring/$newstring/g" test.samples
```
For Squeezemeta the sample file and data should be in same directory 
```
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
```

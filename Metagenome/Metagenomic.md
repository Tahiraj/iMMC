# Taxonomic classification of long reads from shotgun metagenomics sequencing.

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
For single long read use only `pair1`.  
The code below create a directory and new sample file and also replaces the `old sample name` with `new sample name`.

```
#!/bin/bash
#SBATCH --job-name=sequeez
#SBATCH --output=./Log/sequeez.%A_%a.out
#SBATCH --error=./Log/sequeez.%A_%a.err
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --array=1-73

# iMMC/Metagenome/squeezeMeta

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
#SBATCH --mem=360GB
#SBATCH --array=1-73

module load squeezemeta/1.6.3

sample=`cat sample.list | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

cd ${sample}

sqm_longreads.pl -p ${sample} -s test.samples -f data/ -t 32
```

## Extract table for downstream processing

We use [sqmreads2tables.py](https://github.com/jtamames/SqueezeMeta/blob/master/utils/sqmreads2tables.py) to extract table for downstram processing.
```
cat sample.list |head
C1
C2
C3 ...
```
All files in project directory, command look like `sqmreads2tables.py /path/to/project /path/to/output_dir --trusted-functions --force-overwrite`

If file size is large, Do not run in terminal node, rather run as a batch script

```
#!/bin/bash
#SBATCH --job-name=sqmtb
#SBATCH --output=./Log/sqmtb.%J.out
#SBATCH --error=./Log/sqmtb.%J.err
#SBATCH --time=20:00
#SBATCH --nodes=1
#SBATCH --mem=360GB

module load squeezemeta/1.6.3
i=0
while ((i++)); read -r sample
do  
sqmreads2tables.py ${sample}/${sample}/ ${sample}/${sample}_out --trusted-functions --force-overwrite        
done < sample.list
```

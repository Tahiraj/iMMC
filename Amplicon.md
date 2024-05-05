## Taxonomic classification of single reads from  amplicon-targeted 

Create a file `ibex.config`

```
process {

    executor = 'slurm'

    clusterOptions = '--partition=batch --time=12:00:00 --cpus-per-task=32 --mem=50G'

}
```
```
#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --output=./Log/nextflowF.%A.out
#SBATCH --error=./Log/nextflowF.%A.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=1

time date
module load nextflow singularity;
#nextflow run epi2me-labs/wf-metagenomics --fastq data/ -profile singularity -c ibex.config 

nextflow run epi2me-labs/wf-metagenomics --fastq ../data/ --min_len 1200 --min_read_qual 15 --max_len 1800  --database_set "ncbi_16s_18s_28s_ITS" -profile singularity -c ibex.config
```

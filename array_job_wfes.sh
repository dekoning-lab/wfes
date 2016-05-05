#!/usr/bin/env bash

nparams=`wc -l params.txt | cut -f 1 -d ' '`

#SBATCH --job-name WFES
#SBATCH --output=paper_table/wfes_%j.csv
#SBATCH --array=1-22011
#SBATCH --exclude=node[001-004]
#SBATCH -N1

params=`sed $SLURM_ARRAY_TASK_ID'q;d' params.txt`
srun wfes $params

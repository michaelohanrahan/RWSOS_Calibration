#!/bin/bash

#SBATCH -J ORCHESTRA            # name of job
#SBATCH -o orchestra.out            # output file
#SBATCH -N 1              	# total number of nodes requested 
#SBATCH -t 24:00:00	# Run time (hh:mm:ss) - 24 hours
#SBATCH --partition 4vcpu

# activate the conda environment
source $HOME/miniconda3/bin/activate gc-japan
snakemake -s Snakefile_h7 --configfile config/snake_config_singularity.yml --unlock
snakemake -s Snakefile_h7 --directory $PWD --configfile config/snake_config_singularity.yml --cluster "$SNAKE_SUBMIT_JOB" --jobs 5 --cores 1 --latency-wait 30 --wait-for-files --keep-going --rerun-incomplete --forceall --use-conda --conda-prefix /u/winsemi/miniconda3/envs/gc-japan --group-components run_event=10
conda deactivate 

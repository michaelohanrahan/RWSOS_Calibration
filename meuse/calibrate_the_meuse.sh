#!/bin/bash
#SBATCH --job-name=puget_king_cal_2010_18                                               # Job name
#SBATCH --output=/u/dalmijn/code/non_pr/puget/jobs/logging/calib_%j.log                 # Standard output and error log
#SBATCH --time=4-12:00:00                                                               # Job duration (hh:mm:ss)
#SBATCH --partition 48vcpu
#SBATCH --ntasks=48                                                                     # Number of tasks (analyses) to run
#SBATCH --mail-user=brendan.dalmijn@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

# Set variables
export name=calib
export config_file=calib
export rule=all

# Go one directory up to set PWD
cd ..

# Do the thing!
source ~/apps/miniforge3/etc/profile.d/conda.sh
conda activate puget
# Unlock directory
snakemake --unlock -s Snakefile_$name.smk --configfile cfg/$config_file.yml
# Run the workflow!
snakemake $rule --quiet -s Snakefile_$name.smk -c 48 --configfile cfg/$config_file.yml --rerun-incomplete
# Exit

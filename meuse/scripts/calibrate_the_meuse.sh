#!/bin/bash
#SBATCH --job-name=RWSorchestrator                                             # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-logs/h7/calib_%j.log # Standard output and error log
#SBATCH --time=1-12:00:00  # Job duration (hh:mm:ss)
#SBATCH --partition 1vcpu
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

# Set variables
# Go one directory up to set PWD
# cd ..

# Unlock directory
pixi run snakemake --unlock -s "2_Snakefile" --configfile "config/calib.yml"
# Run the workflow!
pixi run snakemake --quiet -s "2_Snakefile" -c 1 --configfile "config/calib.yml" --rerun-incomplete --profile "h7/" --wait-for-files --rerun-triggers mtime

conda deactivate

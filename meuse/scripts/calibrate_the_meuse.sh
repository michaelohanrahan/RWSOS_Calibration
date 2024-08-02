#!/bin/bash
#SBATCH --job-name=RWSorchestrator                                             # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/calib_orch_%j.log # Standard output and error log
#SBATCH --error=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/calib_orch_%j.err  # Standard output and error log
#SBATCH --time=2-00:00:01  # Job duration (hh:mm:ss)
#SBATCH --partition 4pcpu
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

# Set variables
# Go one directory up to set PWD
cd ..
echo "Changed to $(pwd)"
# Unlock directory
pixi run snakemake -s "2_Snakefile.smk" --profile "$HOME/.config/snakemake/slurm/" --unlock
echo "Unlocked directory"
echo "Running snakemake with 4 cores"
# Run the workflow!
pixi run snakemake -s "2_Snakefile.smk" --profile "$HOME/.config/snakemake/slurm/"


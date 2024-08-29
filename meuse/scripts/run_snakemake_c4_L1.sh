#!/bin/bash
#SBATCH --job-name=RWSorchestrator                                           # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/calib_orch_%j.log # Standard output and error log log
#SBATCH --time=15-00:00:00  # Job duration (hh:mm:ss)
#SBATCH --partition 4pcpu
#SBATCH --cpus-per-task=4  # Number of CPUs
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

cd ..
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" -R initial_instate_tomls_L1 --unlock --quiet
pixi run snakemake -s "2_Snakefile.smk" -c 4 --configfile "config/calib.yml" -n -R initial_instate_tomls_L1 --quiet
pixi run snakemake -s "2_Snakefile.smk" -c 4 --profile "$HOME/.config/snakemake/slurm/" -R initial_instate_tomls_L1

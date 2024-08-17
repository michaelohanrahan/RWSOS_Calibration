#!/bin/bash
#SBATCH --job-name=RWSeval                                             # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/eval_%j.log # Standard output and error log log
#SBATCH --time=6:00:00  # Job duration (hh:mm:ss)
#SBATCH --partition 4pcpu
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

cd ..
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --unlock
pixi run snakemake -s "2_Snakefile.smk" -c 1 --configfile "config/calib.yml" -n
pixi run snakemake -s "2_Snakefile.smk" -c 1 --profile "$HOME/.config/snakemake/slurm/" -R evaluate_L0

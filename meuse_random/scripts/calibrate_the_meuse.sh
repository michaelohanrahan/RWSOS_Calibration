#!/bin/bash
#SBATCH --job-name=RWSorc  # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/calib_orch_%j.log # Standard output and error log log
#SBATCH --time=31-12:00:00  # Job duration (hh:mm:ss)
#SBATCH --partition 4pcpu
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

cd ..
echo "Changed to $(pwd)"
pixi run snakemake -s "2_Snakefile.smk" --profile "$HOME/.config/snakemake/slurm/" --unlock
echo "Unlocked directory"
pixi run snakemake -s "2_Snakefile.smk" --profile "$HOME/.config/snakemake/slurm/" --dry-run --quiet
echo "Running snakemake orchestrator"
pixi run snakemake -s "2_Snakefile.smk" -c 4 --profile "$HOME/.config/snakemake/slurm/"


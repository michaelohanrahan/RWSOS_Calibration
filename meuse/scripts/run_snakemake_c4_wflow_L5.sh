#!/bin/bash
#SBATCH --job-name=RWSorchestrator                                       # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/calib_orch_%j.log # Standard output and error log log
#SBATCH --time=3-12:00:00  # Job duration (hh:mm:ss)
#SBATCH --partition 4pcpu
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

cd ..
echo "Changed to $(pwd)"
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --unlock --quiet
echo "cleanup metadata of previous run"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --cleanup-metadata $(cat ./valid_outputs.txt | tr ',' ' ') -n
echo "Unlocked directory"
pixi run snakemake -s "2_Snakefile.smk" --profile "$HOME/.config/snakemake/slurm/" -n -R wflow_L5 --quiet
echo "Running snakemake orchestrator"
pixi run snakemake -s "2_Snakefile.smk" -c 4 --profile "$HOME/.config/snakemake/slurm/" -R wflow_L5

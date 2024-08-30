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
echo "Unlocked directory"
echo "cleanup metadata of previous run"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --cleanup-metadata $(cat ./level5_valid_outputs.txt | tr ',' ' ')
echo "Cleaned up metadata"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --omit-from $(cat ./level5_valid_outputs.txt | tr ',' ' ') -R wflow_L5  --until done_L5 -n
echo "Running snakemake orchestrator"
pixi run snakemake -s "2_Snakefile.smk" -c 4 --profile "slurm/" --omit-from $(cat ./level5_valid_outputs.txt | tr ',' ' ') -R wflow_L5  --until done_L5
